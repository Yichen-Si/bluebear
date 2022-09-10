#include "cramore.h"
#include "utils.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include <math.h>
#include "ibs0.h"
#include "rare_variant_ibs0.h"
#include <iterator>

double gammapdf(float alpha, float beta, float x) {
  return( std::pow(beta, alpha) * std::pow(x,alpha-1) * exp(-beta*x) / tgamma(alpha) );
}
double lgammapdf(float alpha, float beta, float x) {
  return( alpha * log(beta) - lgamma(alpha) + log(x) - x * beta );
}

// goal -- Annotate age estimates. MAP & CME
int32_t AnnotateAge(int32_t argc, char** argv) {
  std::string inVcf, info_f, reg;
  std::string out, outf;
  int32_t verbose = 10000;
  int32_t max_ac = 10, skip = 1;
  int32_t max_age= 1000;
  double cover = 0.95;
  float max_cm = 10.0, min_cm = 0.001;
  std::vector<std::string> tag_list;
  std::string tag_field = "ANN";

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
    LONG_DOUBLE_PARAM("cap",&max_cm,"Maximum ibd length to trust")
    LONG_DOUBLE_PARAM("cup",&min_cm,"Minimum ibd length to trust")
    LONG_MULTI_STRING_PARAM("key-words",&tag_list, "Functional annotation groups to get posterior summary")

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
  outf = out + ".site.age.vcf";

  // Prior
  std::vector<std::vector<double> > prior;
  for (int32_t i = 0; i <= max_ac; ++i) {
    std::vector<double> vec(max_age+1, 0.0);
    prior.push_back(vec);
  }

  // Store the overall distribution of posterior
  // By annotation groups
  std::map<std::string, std::vector<std::vector<double> > > post_dict;
  if (tag_list.size() > 0) {
    for (auto & v : tag_list) {
      std::vector<std::vector<double> > posterior;
      for (int32_t i = 0; i <= max_ac; ++i) {
        std::vector<double> vec(max_age+2, 0.0);
        posterior.push_back(vec);
      }
      post_dict[v] = posterior;
    }
  }

  std::ifstream ifs;
  std::string line;
  std::vector<std::string> words;
  int32_t ac, age, prior_max_age = 0;
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
    if (age > prior_max_age)
      prior_max_age = age;
    if (age > max_age)
      age = max_age;
    prior[ac][age] += p;
  }
  ifs.close();
  notice("Finish reading prior info");
  if (prior_max_age > max_age) {
    max_age = prior_max_age;
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
  // New info headers for age estimates based on averaged IBD/IBS
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "AgeMLE" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=AgeMLE,Number=1,Type=Integer,Description=\"MLE age (only for doubletons)\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
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

  // New info headers for age estimates based on median IBD/IBS
  // TODO need better tag
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "AgeMAP2" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=AgeMAP2,Number=1,Type=Float,Description=\"Posterior mode of the variant age\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "AgeCME2" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=AgeCME2,Number=1,Type=Float,Description=\"Conditional mean estimate of the variant age\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "AgeLowerBound2" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=AgeLowerBound2,Number=1,Type=Float,Description=\"Lower bound of the credible interval of the variant age\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "AgeUpperBound2" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=AgeUpperBound2,Number=1,Type=Float,Description=\"Upper bound of the credible interval of the variant age\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }

  odw.write_hdr();

  notice("Start reading site information from VCF file.");

  float *info_cm = NULL, *med_cm = NULL;
  float *pair_cm = NULL;
  float cm = 0.0, mcm = 0.0;
  int32_t *info_ac = NULL;
  char  *ann = NULL;
  int32_t n_ac = 0, n_cm = 0, n_bt = 0, n_ann = 0;
  int32_t n_md = 0;
  int32_t nVariant = 0;
  float rg = (1.0-cover)/2.0;

  for(int32_t k=0; odr.read(iv); ++k) {

    if (!bcf_is_snp(iv)) {continue;}

    // periodic message to user
    if ( k % verbose == 0 ) {
      notice("Processing %d markers at %s:%d. Recorded %d variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nVariant);
    }
    bool do_avg = true;
    bool do_med = true;
    if (bcf_get_info_float(odr.hdr, iv, "AvgDist_cM", &info_cm, &n_cm) < 0) {do_avg = false;}
    if (bcf_get_info_float(odr.hdr, iv, "MedDist_cM", &med_cm, &n_md) < 0) {do_med = false;}
    if ((!do_avg) && (!do_med)) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (info_ac[0] > max_ac || info_ac[0] < 2) {continue;}

    ac = info_ac[0];
    int32_t have_tag_field = bcf_get_info_string(odr.hdr, iv, tag_field.c_str(), &ann, &n_ann);
    double prob[max_age+1]; // Posterior mean
    double prob_map = 0.0, w_sum = 0.0, p_sum = 0.0, tmp = 0.0;
    float est_map = 0.0, ml = 0.0;
    float lower = -1.0, upper = -1.0;
    int32_t mle = 0;

    if (do_avg) {

      cm = std::min(info_cm[0],max_cm);
      cm = std::max(cm,min_cm);

      for ( uint32_t i = 1; i < prior[ac].size(); ++i ) {
        tmp = gammapdf(2, 0.02*i, cm);
        tmp *= prior[ac][i];
        if (i < prior[ac].size() - 1) {
          if (ml < tmp) {
            ml = tmp;
            mle = i;
          }
          if (prob_map < tmp) {
            prob_map = tmp;
            est_map = (float) i;
          }
        }
        w_sum += i * tmp;
        p_sum += tmp;
        prob[i] = tmp;
      }

      float cme = (float) w_sum / p_sum;
      float n_sum = 0.0;
      for (int32_t i = 1; i <= max_age; ++i) {
        prob[i] /= p_sum;
        n_sum += prob[i];
        if (lower < 0 && n_sum >= rg) {lower = i;}
        if (upper < 0 && n_sum >= 1.0-rg) {upper=i;break;}
      }

      if (have_tag_field >= 0) {
        std::string fun(ann);
        for (auto & tag : tag_list) {
          if ( fun.find(tag) != std::string::npos ) {
            for (int32_t i = 1; i <= max_age; ++i)
              post_dict[tag][ac][i] += prob[i];
            post_dict[tag][ac].back() += 1;
          }
        }
      }

      if (ac == 2) {
        bcf_update_info_int32(odw.hdr, iv, "AgeMLE", &mle, 1);
      }
      bcf_update_info_float(odw.hdr, iv, "AgeMAP", &est_map, 1);
      bcf_update_info_float(odw.hdr, iv, "AgeCME", &cme, 1);
      bcf_update_info_float(odw.hdr, iv, "AgeLowerBound", &lower, 1);
      bcf_update_info_float(odw.hdr, iv, "AgeUpperBound", &upper, 1);
    }

    if (do_med && ac > 2) {

      cm = std::min(med_cm[0],max_cm);
      cm = std::max(cm,min_cm);
      prob_map = 0.0; w_sum = 0.0; p_sum = 0.0; tmp = 0.0;
      est_map = 0.0; ml = 0.0; mle = 0;
      lower = -1.0; upper = -1.0;

      for ( uint32_t i = 1; i < prior[ac].size(); ++i ) {
        tmp = gammapdf(2, 0.02*i, cm);
        tmp *= prior[ac][i];
        if (i < prior[ac].size() - 1) {
          if (ml < tmp) {
            ml = tmp;
            mle = i;
          }
          if (prob_map < tmp) {
            prob_map = tmp;
            est_map = (float) i;
          }
        }
        w_sum += i * tmp;
        p_sum += tmp;
        prob[i] = tmp;
      }

      float cme = (float) w_sum / p_sum;
      float n_sum = 0.0;
      for (int32_t i = 1; i <= max_age; ++i) {
        prob[i] /= p_sum;
        n_sum += prob[i];
        if (lower < 0 && n_sum >= rg) {lower = i;}
        if (upper < 0 && n_sum >= 1.0-rg) {upper=i;break;}
      }

      if (have_tag_field >= 0) {
        std::string fun(ann);
        for (auto & tag : tag_list) {
          if ( fun.find(tag) != std::string::npos ) {
            for (int32_t i = 1; i <= max_age; ++i)
              post_dict[tag][ac][i] += prob[i];
            post_dict[tag][ac].back() += 1;
          }
        }
      }

      bcf_update_info_float(odw.hdr, iv, "AgeMAP2", &est_map, 1);
      bcf_update_info_float(odw.hdr, iv, "AgeCME2", &cme, 1);
      bcf_update_info_float(odw.hdr, iv, "AgeLowerBound2", &lower, 1);
      bcf_update_info_float(odw.hdr, iv, "AgeUpperBound2", &upper, 1);
    }

    // // Composite likelihood
    // int32_t sz = bcf_get_info_float(odr.hdr, iv, "BtwDist", &pair_cm, &n_bt);

    // if (sz > 1) {
    //   prob_map = 0.0; w_sum = 0.0; p_sum = 0.0; tmp = 0.0;
    //   est_map = 0.0;
    //   lower = -1.0; upper = -1.0;
    //   for ( uint32_t i = 1; i < prior[ac].size(); ++i ) {
    //     tmp = 0.0;
    //     for (int32_t it = 0; it < sz; ++it) {
    //       tmp += lgammapdf(2, 0.02*i, pair_cm[it]);
    //     }
    //     tmp = exp(tmp + log(prior[ac][i]));
    //     if (prob_map < tmp) {
    //       prob_map = tmp;
    //       est_map = (float) i;
    //     }
    //     w_sum += i * tmp;
    //     p_sum += tmp;
    //     prob[i] = tmp;
    //   }
    //   cme = (float) w_sum / p_sum;
    //   n_sum = 0.0;
    //   for (int32_t i = 1; i <= max_age; ++i) {
    //     n_sum += prob[i]/p_sum;
    //     if (lower < 0 && n_sum >= rg) {lower = i;}
    //     if (upper < 0 && n_sum >= 1.0-rg) {upper=i;break;}
    //   }
    // }
    // if (sz > 1 || (sz == 1 && ac == 2) ) {
    //   bcf_update_info_float(odw.hdr, iv, "AgeMAP2", &est_map, 1);
    //   bcf_update_info_float(odw.hdr, iv, "AgeCME2", &cme, 1);
    //   bcf_update_info_float(odw.hdr, iv, "AgeLowerBound2", &lower, 1);
    //   bcf_update_info_float(odw.hdr, iv, "AgeUpperBound2", &upper, 1);
    // }

    odw.write(iv);
    nVariant++;
  }

  odw.close();

  // Output the posterior distribution by functional annotation
  // AC Age Tot Prob Tag
  if (tag_list.size() > 0) {
    outf = out + ".post.age.txt";
    FILE * wf;
    wf = fopen(outf.c_str(), "w");
    for (ac = 2; ac <= max_ac; ++ac) {
      for (auto & v : post_dict) {
        if (v.second[ac].back() <= 0.0) {continue;}
        double psum = 0.0;
        for (int32_t a = 1; a <= max_age; ++a) {
          psum += v.second[ac][a];
        }
        if (psum <= 0.0) {continue;}
        for (int32_t a = 1; a <= max_age; ++a) {
          fprintf(wf, "%d\t%d\t%d\t%.2f\t%.4f\t%s\n", ac, (int32_t) v.second[ac].back(), a, v.second[ac][a], v.second[ac][a]/psum, v.first.c_str());
        }
      }
    }
    fclose(wf);
  }

  return 0;
}
