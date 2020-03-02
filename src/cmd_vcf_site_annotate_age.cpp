#include "cramore.h"
#include "utils.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include <math.h>
#include "ibs0.h"
#include "rare_variant_ibs0.h"

double gammapdf(float alpha, float beta, float x) {
  return( std::pow(beta, alpha) * std::pow(x,alpha-1) * exp(-beta*x) / tgamma(alpha) );
}

// goal -- Annotate age estimates. MAP & CME
int32_t AnnotateAge(int32_t argc, char** argv) {
  std::string inVcf, info_f, reg;
  std::string out, outf;
  int32_t verbose = 10000;
  int32_t max_ac = 10, skip = 1;
  int32_t max_age= 1000;
  double cover = 0.95;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
    LONG_STRING_PARAM("info",&info_f, "Input file containing the prior probabilities")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("max-ac",&max_ac,"Maximal allele count to annotate")
    LONG_INT_PARAM("max-age",&max_age,"Maximal age")
    LONG_INT_PARAM("skip",&skip,"lines in info file to skip")
    LONG_INT_PARAM("cover",&cover,"Credible interval")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file prefix")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || info_f.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --info, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // Prior
  std::vector<std::vector<double> > prior;
  for (int32_t i = 0; i <= max_ac; ++i) {
    std::vector<double> vec(max_age+1, 0.0);
    prior.push_back(vec);
  }

  std::ifstream ifs;
  std::string line;
  std::vector<std::string> words;
  int32_t ac, age;
  double p;
  ifs.open(info_f, std::ifstream::in);
  if (!ifs.is_open()) {
    error("Info file cannot be opened");
  }
  if (skip > 0) {
    for (int32_t i = 0; i < skip; ++i) {
      std::getline(ifs, line);
    }
  }
  while(std::getline(ifs, line)) {
    words.clear();
    split(words, "\t", line);
    try {
      ac = std::stoi(words[0]);
      age= std::stoi(words[1]);
      p = std::stod(words[3]);
    }
    catch (const std::invalid_argument& ia) {
      continue;
    }
    if (age > max_age)
      age = max_age;
    prior[ac][age] = p;
  }

  ifs.close();
  notice("Finish reading prior info");

  outf = out + ".site.age.vcf";

  // bcf reader
  std::vector<GenomeInterval> intervals;
  if (!reg.empty())
    parse_intervals(intervals, "", reg);
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();

  // bcf writer (site only)
  BCFOrderedWriter odw(outf.c_str(),0);
  bcf_hdr_t* hnull = bcf_hdr_subset(odr.hdr, 0, 0, 0);
  bcf_hdr_remove(hnull, BCF_HL_FMT, NULL);
  odw.set_hdr(hnull);

  char buffer[65536];
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "AgeMAP" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=AgeMAP,Number=1,Type=Float,Description=\"Posterior mode of the variant age\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "AgeCME" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=AgeCME,Number=1,Type=Float,Description=\"Conditional mean estimate of the variant age\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "AgeLowerBound" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=AgeLowerBound,Number=1,Type=Float,Description=\"Lower bound of the credible interval of the variant age\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "AgeUpperBound" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=AgeUpperBound,Number=1,Type=Float,Description=\"Upper bound of the credible interval of the variant age\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  odw.write_hdr();

  notice("Start reading site information from VCF file.");

  float *info_cm = NULL;
  float cm = 0.0;
  int32_t *info_ac = NULL;
  int32_t n_ac = 0, n_cm = 0;
  int32_t nVariant = 0;
  float rg = (1.0-cover)/2.0;

  for(int32_t k=0; odr.read(iv); ++k) {

    if (!bcf_is_snp(iv)) {continue;}

    // periodic message to user
    if ( k % verbose == 0 ) {
      notice("Processing %d markers at %s:%d. Recorded %d variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nVariant);
    }
    if (bcf_get_info_float(odr.hdr, iv, "AvgDist_cM", &info_cm, &n_cm) < 0) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (info_ac[0] > max_ac || info_ac[0] < 2) {continue;}
    cm = info_cm[0];
    ac = info_ac[0];
    if (cm < 0) {cm = 0.0;} // Never happen

    double prob[max_age+1];
    double prob_map = 0.0, w_sum = 0.0, p_sum = 0.0, tmp = 0.0;
    float est_map = 0.0;
    float lower = -1.0, upper = -1.0;
    for ( uint32_t i = 1; i < prior[ac].size(); ++i ) {
      tmp = gammapdf(2, 0.04*i, cm);
      if (prob_map < tmp) {
        prob_map = tmp;
        est_map = (float) i;
      }
      w_sum += i * tmp;
      p_sum += tmp;
      prob[i] = tmp;
    }
    float cme = (float) w_sum / p_sum;
    float n_sum = 0.0;
    for (uint32_t i = 1; i <= max_age; ++i) {
      n_sum += prob[i]/p_sum;
      if (lower < 0 && n_sum >= rg) {lower = i;}
      if (upper < 0 && n_sum >= 1.0-rg) {upper=i;break;}
    }

    bcf_update_info_float(odw.hdr, iv, "AgeMAP", &est_map, 1);
    bcf_update_info_float(odw.hdr, iv, "AgeCME", &cme, 1);
    bcf_update_info_float(odw.hdr, iv, "AgeLowerBound", &lower, 1);
    bcf_update_info_float(odw.hdr, iv, "AgeUpperBound", &upper, 1);
// std::cout << ac << '\t' << cm << '\t' << est_map << '\t' << cme << '\n';
    odw.write(iv);
    nVariant++;
  }
  odw.close();
  return 0;
}




