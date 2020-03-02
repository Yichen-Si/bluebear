#include "cramore.h"
#include "utils.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include <iomanip>
#include "ibs0.h"
#include "rare_variant_ibs0.h"

// goal -- Annotate PM posterior.
// TODO: general info field update
int32_t AnnotatePM(int32_t argc, char** argv) {
  std::string inVcf, info_f, reg;
  std::string out, outf;
  int32_t verbose = 10000;
  int32_t target_ac = 2, skip = 1;
  double bin = 0.002, maxval = 3.0, upperbound = 2.0;
  int32_t maxb = std::round(maxval/bin);
  int32_t kmersize = 5, offset = 0;
  std::vector<std::string> keywords;
  std::vector<double> th_vec;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
    LONG_STRING_PARAM("info",&info_f, "Input file containing the information to annotate")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("ac",&target_ac,"Allele count to annotate")
    LONG_INT_PARAM("skip",&skip,"lines in info file to skip")
    LONG_INT_PARAM("kmer",&kmersize, "The kmer to be considered. Cannot exceed what has been annotated in the VCF")
    LONG_DOUBLE_PARAM("bin",&bin,"Piece wise constant")
    LONG_MULTI_STRING_PARAM("key-words",&keywords, "Functional annotation to summarize")
    LONG_MULTI_DOUBLE_PARAM("prob-cut",&th_vec, "Threshold to cut pm probability")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file prefix")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  if (th_vec.size() == 0) {
    th_vec.push_back(0.5);
    th_vec.push_back(0.8);
  }

  // Build ancestry info map
  std::map<std::string, std::vector<double> > ref_map;
  std::map<std::string, std::map<std::string, std::vector<double> > > fn_ct;

  std::vector<std::string> type_key{"CpG","TS","TV"};
  std::vector<double> cvec(2+th_vec.size(), 0.0);
  for (auto & v : type_key) {
    std::map<std::string, std::vector<double> > tmp;
    for (auto & u : keywords) {
      tmp[u] = cvec;
    }
    tmp["Tot"] = cvec;
    fn_ct[v] = tmp;
  }

  std::ifstream ifs;
  std::string line;
  std::vector<std::string> words;
  std::string kmer;
  double x, post;
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
      kmer = words[0];
      x = std::stod(words[1]);
      post = std::stod(words[2]);
    }
    catch (const std::invalid_argument& ia) {
      continue;
    }
    int32_t b = (int32_t) std::round(x/bin);
    auto ptr = ref_map.find(kmer);
    if (ptr == ref_map.end()) {
      std::vector<double> vec(maxb+1,0.0);
      for (int32_t i = 0; i <= b; ++i) {
        vec[i] = post;
      }
      ref_map[kmer] = vec;
    } else {
      ptr->second[b] = post;
    }
  }

  ifs.close();
  notice("Finish reading ref info for %d kmer", ref_map.size());

  outf = out;
  if (!reg.empty()) {
    outf += "_"+reg+".site.pm.vcf";
  } else {
    outf += ".site.pm.vcf";
  }

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
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "PMprob" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=PMprob,Number=1,Type=Float,Description=\"Posterior probability of being a parallel mutation\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  odw.write_hdr();

  notice("Start reading site information from VCF file.");

  // mutation types
  std::vector<char> alphabet {'A','T','C','G'};
  std::set<char> alphaset(alphabet.begin(),alphabet.end());
  std::set<std::string> tsmut{"A-G","C-T","G-A","T-C"};
  std::map<char,char> bpair = { {'A','T'}, {'T','A'}, {'C','G'}, {'G','C'}, {'-','-'} };

  float *info_cm = NULL;
  char  *ctx = NULL;
  char  *ann = NULL;
  char  *cpg = NULL;
  float cm = 0.0;
  int32_t *info_ac = NULL;
  int32_t n_ac = 0, n_cm = 0, n_ctx = 0, n_ann = 0, n_cpg = 0;
  int32_t nVariant = 0;

  for(int32_t k=0; odr.read(iv); ++k) {

    if (!bcf_is_snp(iv)) {continue;}

    // periodic message to user
    if ( k % verbose == 0 ) {
      notice("Processing %d markers at %s:%d. Recorded %d variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nVariant);
    }
    if (bcf_get_info_string(odr.hdr, iv, "ANN", &ann, &n_ann) < 0) {continue;}
    if (bcf_get_info_string(odr.hdr, iv, "CpG", &cpg, &n_cpg) < 0) {continue;}
    if (bcf_get_info_string(odr.hdr, iv, "CONTEXT", &ctx, &n_ctx) < 0) {continue;}
    if (bcf_get_info_float(odr.hdr, iv, "AvgDist_cM", &info_cm, &n_cm) < 0) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (info_ac[0] != target_ac) {continue;}
    cm = info_cm[0];
    if (cm < 0) {cm = 0.0;}

    std::string fun(ann);
    std::string type = "TV";
    char ref = iv->d.allele[0][0];
    char alt = iv->d.allele[1][0];
    std::string mtype(1, ref);
    mtype += "-";
    mtype += alt;
    if ( tsmut.find(mtype) != tsmut.end() ) {
      type = "TS";
    }
    if (cpg[0] == 'T') {
      type = "CpG";
    }

    int32_t b = (int32_t) std::round(cm / bin);
    float prob = 0.0;
    if (cm < upperbound) {
      bool flag = 0;
      for (int32_t it = 0; it < kmersize; ++it) {
        if ( alphaset.find(ctx[offset+it]) == alphaset.end() ) {
          flag = 1;
          break;
        }
      }
      if (flag) {continue;}
      std::string kctx(ctx);
      if (kmersize < ((int32_t) kctx.size())) {
        offset = (kctx.size()-kmersize)/2;
        kctx = kctx.substr(offset, kmersize);
      }

      if ( std::string(1, ref) != kctx.substr((kmersize-1)/2, 1) ) {
        notice("Ref. allele != kmer center (%c/%c v.s. %s)", ref, alt, kctx.c_str());
      }
      kctx = kctx.substr(0,(kmersize-1)/2)+mtype+kctx.substr((kmersize-1)/2+1);
      if ( ref == 'T' || ref == 'G' ) {
        std::string fctx = "";
        for (uint32_t it = kctx.size(); it > 0; --it) {
          fctx += bpair[kctx[it-1]];
        }
        kctx = fctx;
      }

      auto ptr = ref_map.find(kctx);
      if (ptr == ref_map.end()) {
        continue;
      }
      prob = ptr->second[b];
    }

    for (auto & v : keywords) {
      if ( fun.find(v) != std::string::npos ) {
        fn_ct[type][v][0] += 1;
        fn_ct[type][v][1] += prob;
        for (uint32_t t = 0; t < th_vec.size(); ++t) {
          if (prob > th_vec[t])
            fn_ct[type][v][2+t] += 1;
        }
      }
    }
    fn_ct[type]["Tot"][0] += 1;
    fn_ct[type]["Tot"][1] += prob;
    for (uint32_t t = 0; t < th_vec.size(); ++t) {
      if (prob > th_vec[t])
        fn_ct[type]["Tot"][2+t] += 1;
    }

    bcf_update_info_float(odw.hdr, iv, "PMprob", &prob, 1);
    odw.write(iv);
    nVariant++;
  }
  odw.close();

  outf = out;
  if (!reg.empty()) {
    outf += "_"+reg+".fn.ct";
  } else {
    outf += ".fn.ct";
  }
  std::ofstream wf(outf.c_str(), std::ios::out);
  for (auto & v : fn_ct) {
    for (auto & u : v.second) {
      wf << target_ac << '\t' << v.first << '\t' << u.first;
      for (auto & x : u.second) {
        wf << '\t' << x;
      }
      wf << '\n';
    }
  }
  wf.close();
  return 0;
}




