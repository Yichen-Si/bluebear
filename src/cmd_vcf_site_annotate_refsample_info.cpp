#include "cramore.h"
#include "utils.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include <iomanip>
#include "ibs0.h"
#include "rare_variant_ibs0.h"

// goal -- For ultra-rare variants in a (small) sample
// annotate AC & carrier ancestry from another reference dataset
int32_t AnnotateRareWithTotalAC(int32_t argc, char** argv) {
  std::string inVcf, inVcf_ref, reg;
  std::string out, outf;
  std::string pop_sample, pop_ref;
  int32_t verbose = 10000;
  int32_t max_rare_ac = 2;
  int32_t ncol1 = 3, ncol2 = 8;
  int32_t use_gt = 1;

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("in-vcf-ref",&inVcf_ref, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
    LONG_STRING_PARAM("pop-info-1",&pop_sample, "Input info file containing sample ancestry")
    LONG_STRING_PARAM("pop-info-2",&pop_ref, "Input info file containing reference individual ancestry")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("max-rare",&max_rare_ac, "Maximal minor allele count to be considered as anchor for IBD")
    LONG_INT_PARAM("use-gt",&use_gt, "If GT fields are to be used for ref. VCF/BCF")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || pop_sample.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --out, --pop-info-1 are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }
  if (use_gt == 1 && (inVcf_ref.empty()||pop_ref.empty()) ) {
    error("[E:%s:%d %s] --in-vcf-ref and --pop-info-2 are required if use-gt = 1",__FILE__,__LINE__,__FUNCTION__);
  }

  // Build ancestry info map
  std::vector<std::string> super_pop{"AFR","AMR","EAS","EUR","SAS"};
  std::map<std::string, std::string> popinfo;
  std::ifstream ifs;
  std::string line;
  std::vector<std::string> words;
  std::string id, pop;

  if (use_gt == 1) {
    ifs.open(pop_ref, std::ifstream::in);
    if (!ifs.is_open()) {
      error("Ancestry info for reference individuals cannot be opened");
    }
    while(std::getline(ifs, line)) {
      words.clear();
      split(words, "\t", line);
      try {
        id = words[0];
        pop= words[ncol2-1];
      }
      catch (const std::invalid_argument& ia) {
        continue;
      }
      popinfo[id] = pop;
    }
    ifs.close();
  }

  ifs.open(pop_sample, std::ifstream::in);
  if (!ifs.is_open()) {
    error("Ancestry info for sample individuals cannot be opened");
  }
  while(std::getline(ifs, line)) {
    words.clear();
    split(words, "\t", line);
    try {
      id = words[0];
      pop= words[ncol1-1];
    }
    catch (const std::invalid_argument& ia) {
      continue;
    }
    popinfo[id] = pop;
  }
  ifs.close();

  outf = out;
  if (!reg.empty()) {
    outf += "_"+reg+".site.acinfo.vcf";
  } else {
    outf += ".site.acinfo.vcf";
  }

  // bcf reader
  std::vector<GenomeInterval> intervals;
  if (!reg.empty())
    parse_intervals(intervals, "", reg);
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();
  BCFOrderedReader odr2(inVcf_ref, intervals);
  bcf1_t* iv2 = bcf_init();


  // bcf writer (site only)
  BCFOrderedWriter odw(outf.c_str(),0);
  bcf_hdr_t* hnull = bcf_hdr_subset(odr.hdr, 0, 0, 0);
  bcf_hdr_remove(hnull, BCF_HL_FMT, NULL);
  odw.set_hdr(hnull);
  char buffer[65536];
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "CarrierPop" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=CarrierPop,Number=1,Type=String,Description=\"Super population of the carriers\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "Carrier" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=Carrier,Number=1,Type=String,Description=\"ID of the carriers\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "TotAC" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=TotAC,Number=2,Type=Integer,Description=\"AC,AN in a larger sample\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if (use_gt == 1) {
    if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "SharedPop" ) < 0 ) {
      sprintf(buffer,"##INFO=<ID=SharedPop,Number=1,Type=String,Description=\"Super population of the carriers in the larger sample\">\n");
      bcf_hdr_append(odw.hdr, buffer);
    }
    if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "MatchedPop" ) < 0 ) {
      sprintf(buffer,"##INFO=<ID=MatchedPop,Number=2,Type=Integer,Description=\"x,y: x more carriers from the same super pop and y from others (only for singletons) \">\n");
      bcf_hdr_append(odw.hdr, buffer);
    }
  }
  odw.write_hdr();

  // handle filter string
  std::string filter_str;
  int32_t filter_logic = 0;
  if ( vfilt.include_expr.empty() ) {
    if ( vfilt.exclude_expr.empty() ) {
      // do nothing
    }
    else {
      filter_str = vfilt.exclude_expr;
      filter_logic |= FLT_EXCLUDE;
    }
  }
  else {
    if ( vfilt.exclude_expr.empty() ) {
      filter_str = vfilt.include_expr;
      filter_logic |= FLT_INCLUDE;
    }
    else {
      error("[E:%s:%d %s] Cannot use both --include-expr and --exclude-expr options",__FILE__,__LINE__,__FUNCTION__);
    }
  }

  filter_t* filt = NULL;
  if ( filter_logic != 0 )
    filter_init(odr.hdr, filter_str.c_str());

  // handle --apply-filtrs
  std::vector<int32_t> req_flt_ids;
  if ( !vfilt.required_filters.empty() ) {
    for(int32_t i=0; i < (int32_t)vfilt.required_filters.size(); ++i) {
      req_flt_ids.push_back(bcf_hdr_id2int(odr.hdr, BCF_DT_ID, vfilt.required_filters[i].c_str()));
    }
  }

  int32_t nVariant = 0;
  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
  int32_t nsamples2 = bcf_hdr_nsamples(odr2.hdr);

  notice("Started Reading site information from VCF file, identifying %d samples", nsamples);
  if ( nsamples == 0 )
    error("FATAL ERROR: The VCF does not have any samples with genotypes");

  int32_t *p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t *info_ac = NULL;
  int32_t *info_an = NULL;
  int32_t n_ac = 0, n_an = 0, pos = 0;
  odr2.read(iv2);
  bcf_unpack(iv2, BCF_UN_FLT);
  bool flag = 0;

  for(int32_t k=0; odr.read(iv); ++k) {

    if (!bcf_is_snp(iv)) {continue;}

    if (iv->pos+1 == pos) {
      continue; // Jump over tri-allelic sites
    }
    // periodic message to user
    if ( k % verbose == 0 )
      notice("Processing %d markers at %s:%d. Recorded %d variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nVariant);

    // unpack FILTER column
    bcf_unpack(iv, BCF_UN_FLT);
    // check --apply-filters
    bool has_filter = req_flt_ids.empty() ? true : false;
    if ( ! has_filter ) {
      for(int32_t i=0; i < iv->d.n_flt; ++i) {
        for(int32_t j=0; j < (int32_t)req_flt_ids.size(); ++j) {
          if ( req_flt_ids[j] == iv->d.flt[i] )
            has_filter = true;
        }
      }
    }
    if ( ! has_filter ) { continue; }
    // check filter logic
    if ( filt != NULL ) {
      int32_t ret = filter_test(filt, iv, NULL);
      if ( filter_logic == FLT_INCLUDE ) { if ( !ret)  has_filter = false; }
      else if ( ret ) { has_filter = false; }
    }
    if ( ! has_filter ) { continue; }

    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
    int32_t ac = info_ac[0], an = info_an[0];
    bool flip = 0;
    if (ac > max_rare_ac && an - ac > max_rare_ac) {continue;}
    if (ac == 0 || an == ac) {continue;}
    if (ac > an/2) {
      ac = an - ac;
      flip = 1;
    }
    int32_t org_ac = ac;
    std::string org_pop;
    std::string org_ref = iv->d.allele[0];
    std::string org_alt = iv->d.allele[1];

    // extract genotype and apply genotype level filter
    if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }

    std::stringstream ss;
    std::stringstream ss2;
    int32_t ct = 0;
    for(int32_t i=0; i < nsamples; ++i) {
      int32_t g1 = p_gt[2*i];
      int32_t g2 = p_gt[2*i+1];
      int32_t geno;
      if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
        continue;
      } else {
        geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
        if (flip) {
          geno = ((bcf_gt_allele(g1) == 0) ? 1 : 0) + ((bcf_gt_allele(g2) == 0) ? 1 : 0);
        }
        if (geno >= 1) {
          std::string tmp(odr.hdr->samples[i]);
          if (ct > 0) {ss2 << ",";}
          ss2 << tmp;
          auto ptr = popinfo.find(tmp);
          if (ptr == popinfo.end()) {
            continue;
          }
          tmp = ptr->second;
          if (ct > 0) {ss << ",";}
          else {org_pop = tmp;}
          ss << tmp;
          ct++;
        }
      }
    } // Found carriers
    if (ct == 0 || org_pop.empty()) {
      continue;
    }

    bcf1_t* nv = bcf_dup(iv);
    bcf_unpack(nv, BCF_UN_FLT);
    bcf_subset(odw.hdr, nv, 0, 0);

    bcf_update_info_string(odw.hdr, nv, "CarrierPop", ss.str().c_str());
    bcf_update_info_string(odw.hdr, nv, "Carrier", ss2.str().c_str());

    // Find corresponding position in ref
    if (iv2->pos > iv->pos) {
      continue;
    }
    while(iv2->pos < iv->pos) {
      if (!odr2.read(iv2)) {
        flag = 1;
        break;
      }
      bcf_unpack(iv2, BCF_UN_FLT);
    }
    if (flag) {break;}
    if (iv2->pos >  iv->pos) { // ref does not contain this variant
      // never happen if ref contains the smaller sample
      int32_t totac[2] = {0,2*nsamples2};
      bcf_update_info_int32(odw.hdr, nv, "TotAC", totac, 2);
      continue;
    }
    while (iv2->pos == iv->pos) {
      bcf_unpack(iv2, BCF_UN_FLT);
      std::string tmp_ref = iv2->d.allele[0];
      std::string tmp_alt = iv2->d.allele[1];
      if (tmp_ref == org_ref && tmp_alt == org_alt) {
        break;
      } else {
        if (!odr2.read(iv2)) {
          flag = 1;
          break;
        }
      }
    }
    if (flag) {break;}

    if (iv2->pos > iv->pos) {
      int32_t totac[2] = {0,2*nsamples2};
      bcf_update_info_int32(odw.hdr, nv, "TotAC", totac, 2);
      continue;
    }

    if (iv2->pos == iv->pos) {

      if (bcf_get_info_int32(odr2.hdr, iv2, "AC", &info_ac, &n_ac) < 0) {continue;}
      if (bcf_get_info_int32(odr2.hdr, iv2, "AN", &info_an, &n_an) < 0) {continue;}
      ac = info_ac[0], an = info_an[0];
      flip = 0;
      if (ac > an/2) {
        ac = an - ac;
        flip = 1;
      }

      int32_t totac[2] = {ac,an};
      bcf_update_info_int32(odw.hdr, nv, "TotAC", totac, 2);

if (use_gt == 1) {
      // extract genotype and apply genotype level filter
      if ( bcf_get_genotypes(odr2.hdr, iv2, &p_gt, &n_gt) < 0 ) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr2.hdr, iv2->rid), iv2->pos+1);
      }

      std::map<std::string, int32_t> refcarry;
      std::map<std::string, int32_t> refadmix;
      for (auto & v : super_pop) {
        refcarry[v] = 0;
      }

      std::stringstream().swap(ss);

      ct = 0;
      for(int32_t i=0; i < nsamples2; ++i) {
        int32_t g1 = p_gt[2*i];
        int32_t g2 = p_gt[2*i+1];
        int32_t geno;
        if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
          continue;
        } else {
          geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
          if (flip) {
            geno = ((bcf_gt_allele(g1) == 0) ? 1 : 0) + ((bcf_gt_allele(g2) == 0) ? 1 : 0);
          }
          if (geno >= 1) {
            std::string tmp(odr2.hdr->samples[i]);
            auto ptr = popinfo.find(tmp);
            if (ptr == popinfo.end()) {
              refadmix["NA"]++;
              continue;
            }
            tmp = ptr->second;
            auto ptr2 = refcarry.find(tmp);
            if (ptr2 == refcarry.end()) {
              refadmix[tmp]++;
            } else {
              ptr2->second ++;
            }
            ct++;
          }
        }
      } // Found carriers

      for (auto & v : super_pop) {
        ss << v << ":" << refcarry[v] << ",";
      }
      if (refadmix.size() > 0) {
        for (auto & it : refadmix) {
          ss << "(" << it.first << "):" << it.second << ",";
        }
      }
      std::string sstr = ss.str(); sstr.pop_back();
      bcf_update_info_string(odw.hdr, nv, "SharedPop", sstr.c_str());

      if (org_ac == 1) {
        int32_t x = 0, y = 0;
        auto ptr = refcarry.find(org_pop);
        if (ptr != refcarry.end()) {
          x = ptr->second;
        }
        y = ac - x;
        int32_t nmatch[2] = {x,y};
        bcf_update_info_int32(odw.hdr, nv, "MatchedPop", nmatch, 2);
      }
}


    }

    nVariant++;
    odw.write(nv);
    bcf_destroy(nv);
  }
  odw.close();
  return 0;
}














// =========================================================
// (Downstream summary of the output of the above function)
// =========================================================

// goal -- Summarize the AC in a large ref. pop of singletons in a smaller sample
int32_t cmdVcfSingletonRefAC(int32_t argc, char** argv) {

  std::string inVcf, reg;
  std::string out, outf;
  int32_t verbose = 10000;
  int32_t full_pop = 0;

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("full-pop",&full_pop, "If reference carriers' ancestry is to be considered")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // Build ancestry info map
  std::vector<std::string> super_pop{"AFR","AMR","EAS","EUR","SAS"};
  // int32_t npop = (int32_t) super_pop.size();
  std::map<std::string, std::array<int32_t, 2> > singlecount;
  std::map<std::string, std::map<int32_t, int32_t> * > totcount;
  std::map<std::string, std::map<int32_t, int32_t> * > popcount;
  for ( auto & v : super_pop ) {
    singlecount[v][0] = 0; singlecount[v][1] = 0;
    totcount[v] = new std::map<int32_t, int32_t>;
    popcount[v] = new std::map<int32_t, int32_t>;
  }

  // bcf reader
  std::vector<GenomeInterval> intervals;
  if (!reg.empty())
    parse_intervals(intervals, "", reg);
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();

  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "CarrierPop" ) < 0 || bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "TotAC" ) < 0 ) {
    error("Cannot find required info field CarrierPop and TotAC");
  }
  if (full_pop == 1 && bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "MatchedPop" ) < 0) {
    error("Cannot find info field MatchedPop");
  }

  // handle filter string
  std::string filter_str;
  int32_t filter_logic = 0;
  if ( vfilt.include_expr.empty() ) {
    if ( vfilt.exclude_expr.empty() ) {
      // do nothing
    }
    else {
      filter_str = vfilt.exclude_expr;
      filter_logic |= FLT_EXCLUDE;
    }
  }
  else {
    if ( vfilt.exclude_expr.empty() ) {
      filter_str = vfilt.include_expr;
      filter_logic |= FLT_INCLUDE;
    }
    else {
      error("[E:%s:%d %s] Cannot use both --include-expr and --exclude-expr options",__FILE__,__LINE__,__FUNCTION__);
    }
  }

  filter_t* filt = NULL;
  if ( filter_logic != 0 )
    filter_init(odr.hdr, filter_str.c_str());

  // handle --apply-filtrs
  std::vector<int32_t> req_flt_ids;
  if ( !vfilt.required_filters.empty() ) {
    for(int32_t i=0; i < (int32_t)vfilt.required_filters.size(); ++i) {
      req_flt_ids.push_back(bcf_hdr_id2int(odr.hdr, BCF_DT_ID, vfilt.required_filters[i].c_str()));
    }
  }

  int32_t nVariant = 0, pos = 0;
  int32_t *info_ac = NULL, *tot_ac = NULL;
  int32_t *m_ac = NULL;
  int32_t n_ac = 0, n_tac = 0, n_mac = 0, n_pop = 0;
  char *org_pop = NULL;

  for(int32_t k=0; odr.read(iv); ++k) {

    if (!bcf_is_snp(iv)) {continue;}
    if (iv->pos == pos) {continue;}
    pos = iv->pos;
    // periodic message to user
    if ( k % verbose == 0 )
      notice("Processing %d markers at %s:%d. Recorded %d singletons", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nVariant);

    // unpack FILTER column
    bcf_unpack(iv, BCF_UN_FLT);
    // check --apply-filters
    bool has_filter = req_flt_ids.empty() ? true : false;
    if ( ! has_filter ) {
      for(int32_t i=0; i < iv->d.n_flt; ++i) {
        for(int32_t j=0; j < (int32_t)req_flt_ids.size(); ++j) {
          if ( req_flt_ids[j] == iv->d.flt[i] )
            has_filter = true;
        }
      }
    }
    if ( ! has_filter ) { continue; }
    // check filter logic
    if ( filt != NULL ) {
      int32_t ret = filter_test(filt, iv, NULL);
      if ( filter_logic == FLT_INCLUDE ) { if ( !ret)  has_filter = false; }
      else if ( ret ) { has_filter = false; }
    }
    if ( ! has_filter ) { continue; }

    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    int32_t ac = info_ac[0];
    if (ac != 1) {continue;}

    if (bcf_get_info_string(odr.hdr, iv, "CarrierPop", &org_pop, &n_pop) < 0) {continue;}
    std::string opop(org_pop);

    if (bcf_get_info_int32(odr.hdr, iv, "TotAC", &tot_ac, &n_tac) < 0) {continue;}
    int32_t t_ac = tot_ac[0];
    if (t_ac == 1) {
      singlecount[opop][0]++;
    }
    singlecount[opop][1]++;
    totcount[opop]->operator[](t_ac)++;

    if (full_pop == 1) {
      if (bcf_get_info_int32(odr.hdr, iv, "MatchedPop", &m_ac, &n_mac) < 0) {continue;}
      int32_t m = m_ac[0];
      popcount[opop]->operator[](m)++;
    }

    nVariant++;
  }

  if (!reg.empty()) {
    out += "_"+reg;
  }
  FILE* wf;

  outf = out + "_tot.count";
  wf = fopen(outf.c_str(), "w");
  for (auto & it : totcount) {
    for (auto & v : (*it.second)) {
      fprintf(wf, "%s\t%d\t%d\tTot\n", it.first.c_str(), v.first, v.second);
    }
  }
  if (full_pop == 1) {
    for (auto & it : popcount) {
      for (auto & v : (*it.second)) {
        fprintf(wf, "%s\t%d\t%d\tSamePop\n", it.first.c_str(), v.first, v.second);
      }
    }
  }
  fclose(wf);

  outf = out + "_single.count";
  wf = fopen(outf.c_str(), "w");
  for (auto & v : singlecount) {
    fprintf(wf, "%s\t%d\t%d\n", v.first.c_str(), v.second[0], v.second[1]);
  }
  fclose(wf);

  return 0;
}









// =========================================================
// A more general function. Test.
// =========================================================

// goal -- For (ultra-)rare variants in a (small) sample
// annotate AC from another reference dataset
int32_t AnnotateRefAC(int32_t argc, char** argv) {
  std::string inVcf, inVcf_ref, reg;
  std::string out, outf;
  std::string pop_sample;
  std::string pop_ref;
  int32_t verbose = 10000;
  int32_t max_ac = 10;
  double  max_af = -1.0;
  int32_t ncol1 = 3, ncol2 = 8;
  int32_t superpop = 1, carry = 0;
  int32_t use_gt = 0, keepall = 0;

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("in-vcf-ref",&inVcf_ref, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
    LONG_STRING_PARAM("pop-info-1",&pop_sample, "Input info file containing sample ancestry")
    LONG_STRING_PARAM("pop-info-2",&pop_ref, "Input info file containing reference individual ancestry")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("max-ac",&max_ac, "Maximal minor allele count to be annotated")
    LONG_DOUBLE_PARAM("max-af",&max_af, "Maximal minor allele frequency to be annotated")
    LONG_INT_PARAM("use-gt",&use_gt, "If GT fields are to be used for ref. VCF/BCF")
    LONG_INT_PARAM("carrier-id",&carry, "If carriers' id is to be annotated")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")
    LONG_INT_PARAM("keep-all",&keepall, "If variants (regardless of AC) are to be output")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }
  if (use_gt == 1 && (inVcf_ref.empty()||pop_ref.empty()) ) {
    error("[E:%s:%d %s] --in-vcf-ref and --pop-info-2 are required if use-gt = 1",__FILE__,__LINE__,__FUNCTION__);
  }
  if (pop_sample.empty())
    superpop = 0;

  // Build ancestry info map
  std::vector<std::string> super_pop{"AFR","AMR","EAS","EUR","SAS"};
  std::map<std::string, std::string> popinfo;
  std::ifstream ifs;
  std::string line;
  std::vector<std::string> words;
  std::string id, pop;

  if (use_gt == 1) {
    ifs.open(pop_ref, std::ifstream::in);
    if (!ifs.is_open()) {
      error("Ancestry info for reference individuals cannot be opened");
    }
    while(std::getline(ifs, line)) {
      words.clear();
      split(words, "\t", line);
      try {
        id = words[0];
        pop= words[ncol2-1];
      }
      catch (const std::invalid_argument& ia) {
        continue;
      }
      popinfo[id] = pop;
    }
    ifs.close();
  }

  if (superpop == 1) {
    ifs.open(pop_sample, std::ifstream::in);
    if (!ifs.is_open()) {
      error("Ancestry info for sample individuals cannot be opened");
    }
    while(std::getline(ifs, line)) {
      words.clear();
      split(words, "\t", line);
      try {
        id = words[0];
        pop= words[ncol1-1];
      }
      catch (const std::invalid_argument& ia) {
        continue;
      }
      popinfo[id] = pop;
    }
    ifs.close();
  }

  outf = out;
  if (!reg.empty()) {
    outf += "_"+reg+".site.refac.vcf";
  } else {
    outf += ".site.refac.vcf";
  }

  // bcf reader
  std::vector<GenomeInterval> intervals;
  if (!reg.empty())
    parse_intervals(intervals, "", reg);
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();
  BCFOrderedReader odr2(inVcf_ref, intervals);
  bcf1_t* iv2 = bcf_init();

  // bcf writer (site only)
  BCFOrderedWriter odw(outf.c_str(),0);
  bcf_hdr_t* hnull = bcf_hdr_subset(odr.hdr, 0, 0, 0);
  bcf_hdr_remove(hnull, BCF_HL_FMT, NULL);
  odw.set_hdr(hnull);

  char buffer[65536];
  if ( superpop == 1 && bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "Pop" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=Pop,Number=1,Type=String,Description=\"Super population of the carriers\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( carry == 1 && bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "Carrier" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=Carrier,Number=1,Type=String,Description=\"ID of the carriers\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "TotAC" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=TotAC,Number=2,Type=Integer,Description=\"AC,AN in a larger sample\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if (use_gt == 1) {
    if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "SharedPop" ) < 0 ) {
      sprintf(buffer,"##INFO=<ID=SharedPop,Number=1,Type=String,Description=\"Super population of the carriers in the larger sample\">\n");
      bcf_hdr_append(odw.hdr, buffer);
    }
    if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "MatchedPop" ) < 0 ) {
      sprintf(buffer,"##INFO=<ID=MatchedPop,Number=2,Type=Integer,Description=\"x,y: x more carriers from the same super pop and y from others (only for singletons) \">\n");
      bcf_hdr_append(odw.hdr, buffer);
    }
  }
  odw.write_hdr();

  // handle filter string
  std::string filter_str;
  int32_t filter_logic = 0;
  if ( vfilt.include_expr.empty() ) {
    if ( vfilt.exclude_expr.empty() ) {
      // do nothing
    }
    else {
      filter_str = vfilt.exclude_expr;
      filter_logic |= FLT_EXCLUDE;
    }
  }
  else {
    if ( vfilt.exclude_expr.empty() ) {
      filter_str = vfilt.include_expr;
      filter_logic |= FLT_INCLUDE;
    }
    else {
      error("[E:%s:%d %s] Cannot use both --include-expr and --exclude-expr options",__FILE__,__LINE__,__FUNCTION__);
    }
  }

  filter_t* filt = NULL;
  if ( filter_logic != 0 )
    filter_init(odr.hdr, filter_str.c_str());

  // handle --apply-filtrs
  std::vector<int32_t> req_flt_ids;
  if ( !vfilt.required_filters.empty() ) {
    for(int32_t i=0; i < (int32_t)vfilt.required_filters.size(); ++i) {
      req_flt_ids.push_back(bcf_hdr_id2int(odr.hdr, BCF_DT_ID, vfilt.required_filters[i].c_str()));
    }
  }

  int32_t nVariant = 0;
  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
  int32_t nsamples2 = bcf_hdr_nsamples(odr2.hdr);

  notice("Started Reading site information from VCF file, identifying %d samples", nsamples);
  if ( nsamples == 0 ) {
    carry = 0; superpop = 0;
    notice("The VCF does not have any samples with genotypes");
  }
  if ( nsamples2 == 0 ) {use_gt = 0;}
  if (max_af > 0.0 && nsamples > 0) {
    max_ac = std::max(max_ac, (int32_t) (2*nsamples*max_af));
  }

  int32_t *p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t *info_ac = NULL;
  int32_t *info_an = NULL;
  int32_t n_ac = 0, n_an = 0;
  odr2.read(iv2);
  bcf_unpack(iv2, BCF_UN_FLT);
  bool flag = 0;

  for(int32_t k=0; odr.read(iv); ++k) {

    if (!bcf_is_snp(iv)) {continue;}

    // periodic message to user
    if ( k % verbose == 0 )
      notice("Processing %d markers at %s:%d. Recorded %d variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nVariant);

    // unpack FILTER column
    bcf_unpack(iv, BCF_UN_FLT);
    // check --apply-filters
    bool has_filter = req_flt_ids.empty() ? true : false;
    if ( ! has_filter ) {
      for(int32_t i=0; i < iv->d.n_flt; ++i) {
        for(int32_t j=0; j < (int32_t)req_flt_ids.size(); ++j) {
          if ( req_flt_ids[j] == iv->d.flt[i] )
            has_filter = true;
        }
      }
    }
    if ( ! has_filter ) { continue; }
    // check filter logic
    if ( filt != NULL ) {
      int32_t ret = filter_test(filt, iv, NULL);
      if ( filter_logic == FLT_INCLUDE ) { if ( !ret)  has_filter = false; }
      else if ( ret ) { has_filter = false; }
    }
    if ( ! has_filter ) { continue; }

    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
    int32_t ac = info_ac[0], an = info_an[0];
    bool flip = 0;
    if (ac > max_ac && an - ac > max_ac) {
      if (keepall) {
        bcf1_t* nv = bcf_dup(iv);
        bcf_unpack(nv, BCF_UN_FLT);
        bcf_subset(odw.hdr, nv, 0, 0);
        odw.write(nv);
        bcf_destroy(nv);
      }
      continue;
    }
    if (ac == 0 || an == ac) {continue;}
    if (ac > an/2) {
      ac = an - ac;
      flip = 1;
    }
    int32_t ct = 0;
    std::string org_pop;
    std::string org_ref = iv->d.allele[0];
    std::string org_alt = iv->d.allele[1];
    std::stringstream ss;
    std::string sstr;

    bcf1_t* nv = bcf_dup(iv);
    bcf_unpack(nv, BCF_UN_FLT);
    bcf_subset(odw.hdr, nv, 0, 0);

    if (superpop == 1) {
      // extract genotype and apply genotype level filter
      if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
      }

      std::stringstream ss2;
      std::map<std::string, int32_t> popct;
      for (auto & v : super_pop) {
        popct[v] = 0;
      }

      ct = 0;
      for(int32_t i=0; i < nsamples; ++i) {
        int32_t g1 = p_gt[2*i];
        int32_t g2 = p_gt[2*i+1];
        int32_t geno;
        if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
          continue;
        } else {
          geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
          if (flip) {
            geno = ((bcf_gt_allele(g1) == 0) ? 1 : 0) + ((bcf_gt_allele(g2) == 0) ? 1 : 0);
          }
          if (geno >= 1) {
            std::string tmp(odr.hdr->samples[i]);
            auto ptr = popinfo.find(tmp);
            if (ptr != popinfo.end()) {
              popct[ptr->second]++;
              if (ac == 1 && ct == 0)
                org_pop = ptr->second;
            }
            if (ct > 0) {ss2 << ",";}
            ss2 << i;
            ct++;
          }
        }
      } // Found carriers

      for (auto & v : super_pop) {
        ss << v << ":" << popct[v] << ",";
      }
      sstr = ss.str(); sstr.pop_back();
      bcf_update_info_string(odw.hdr, nv, "Pop", sstr.c_str());
      if (carry == 1)
        bcf_update_info_string(odw.hdr, nv, "Carrier", ss2.str().c_str());
    }

    // Find corresponding position in ref
    while(iv2->pos < iv->pos) {
      if (!odr2.read(iv2)) {
        flag = 1;
        break;
      }
      bcf_unpack(iv2, BCF_UN_FLT);
    }
    if (flag) {break;}
    while (iv2->pos == iv->pos) {
      bcf_unpack(iv2, BCF_UN_FLT);
      std::string tmp_ref = iv2->d.allele[0];
      std::string tmp_alt = iv2->d.allele[1];
      if (tmp_ref == org_ref && tmp_alt == org_alt) {
        break;
      } else {
        if (!odr2.read(iv2)) {
          flag = 1;
          break;
        }
      }
    }
    if (flag) {break;}
    if (iv2->pos > iv->pos) { // ref does not contain this variant
      // never happen if ref contains the smaller sample
      int32_t totac[2] = {0,0};
      bcf_update_info_int32(odw.hdr, nv, "TotAC", totac, 2);
    }
    else if (iv2->pos == iv->pos) {

      if (bcf_get_info_int32(odr2.hdr, iv2, "AC", &info_ac, &n_ac) < 0) {continue;}
      if (bcf_get_info_int32(odr2.hdr, iv2, "AN", &info_an, &n_an) < 0) {continue;}
      ac = info_ac[0], an = info_an[0];
      flip = 0;
      if (ac > an/2) {
        ac = an - ac;
        flip = 1;
      }

      int32_t totac[2] = {ac,an};
      bcf_update_info_int32(odw.hdr, nv, "TotAC", totac, 2);

      if (use_gt == 1) {
        // extract genotype and apply genotype level filter
        if ( bcf_get_genotypes(odr2.hdr, iv2, &p_gt, &n_gt) < 0 ) {
          error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr2.hdr, iv2->rid), iv2->pos+1);
        }

        std::map<std::string, int32_t> refcarry;
        std::map<std::string, int32_t> refadmix;
        for (auto & v : super_pop) {
          refcarry[v] = 0;
        }
        std::stringstream().swap(ss);
        ct = 0;
        for(int32_t i=0; i < nsamples2; ++i) {
          int32_t g1 = p_gt[2*i];
          int32_t g2 = p_gt[2*i+1];
          int32_t geno;
          if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
            continue;
          } else {
            geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
            if (flip) {
              geno = ((bcf_gt_allele(g1) == 0) ? 1 : 0) + ((bcf_gt_allele(g2) == 0) ? 1 : 0);
            }
            if (geno >= 1) {
              std::string tmp(odr2.hdr->samples[i]);
              auto ptr = popinfo.find(tmp);
              if (ptr == popinfo.end()) {
                refadmix["NA"]++;
                continue;
              }
              tmp = ptr->second;
              auto ptr2 = refcarry.find(tmp);
              if (ptr2 == refcarry.end()) {
                refadmix[tmp]++;
              } else {
                ptr2->second ++;
              }
              ct++;
            }
          }
        } // Found carriers

        for (auto & v : super_pop) {
          ss << v << ":" << refcarry[v] << ",";
        }
        if (refadmix.size() > 0) {
          for (auto & it : refadmix) {
            ss << "(" << it.first << "):" << it.second << ",";
          }
        }
        sstr = ss.str(); sstr.pop_back();
        bcf_update_info_string(odw.hdr, nv, "SharedPop", sstr.c_str());

        if (!org_pop.empty()) {
          int32_t x = 0, y = 0;
          auto ptr = refcarry.find(org_pop);
          if (ptr != refcarry.end()) {
            x = ptr->second;
          }
          y = ac - x;
          int32_t nmatch[2] = {x,y};
          bcf_update_info_int32(odw.hdr, nv, "MatchedPop", nmatch, 2);
        }
      }
    } else { // Should never happen
    }

    nVariant++;
    odw.write(nv);
    bcf_destroy(nv);
  }
  odw.close();
  return 0;
}







// =========================================================
// (Downstream summary of the output of the above function)
// =========================================================

// goal -- Summarize the AC in a large ref. pop of singletons in a smaller sample
int32_t cmdVcfSummerRefAC(int32_t argc, char** argv) {

  std::string inVcf, reg;
  std::string out, outf;
  int32_t verbose = 10000;
  int32_t max_ac  = 10;

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("max-ac",&max_ac, "Maximal minor allele count to be included")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // Build ancestry info map
  std::vector<std::string> super_pop{"AFR","AMR","EAS","EUR","SAS"};
  // int32_t npop = (int32_t) super_pop.size();
  std::map<std::string, int32_t* > totcount;
  for ( auto & v : super_pop ) {
    totcount[v] = new int32_t[(max_ac+1)*2]{0};
  }

  // bcf reader
  std::vector<GenomeInterval> intervals;
  if (!reg.empty())
    parse_intervals(intervals, "", reg);
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();

  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "Pop" ) < 0 || bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "TotAC" ) < 0 ) {
    error("Cannot find required info field Pop and TotAC");
  }

  // handle filter string
  std::string filter_str;
  int32_t filter_logic = 0;
  if ( vfilt.include_expr.empty() ) {
    if ( vfilt.exclude_expr.empty() ) {
      // do nothing
    }
    else {
      filter_str = vfilt.exclude_expr;
      filter_logic |= FLT_EXCLUDE;
    }
  }
  else {
    if ( vfilt.exclude_expr.empty() ) {
      filter_str = vfilt.include_expr;
      filter_logic |= FLT_INCLUDE;
    }
    else {
      error("[E:%s:%d %s] Cannot use both --include-expr and --exclude-expr options",__FILE__,__LINE__,__FUNCTION__);
    }
  }

  filter_t* filt = NULL;
  if ( filter_logic != 0 )
    filter_init(odr.hdr, filter_str.c_str());

  // handle --apply-filtrs
  std::vector<int32_t> req_flt_ids;
  if ( !vfilt.required_filters.empty() ) {
    for(int32_t i=0; i < (int32_t)vfilt.required_filters.size(); ++i) {
      req_flt_ids.push_back(bcf_hdr_id2int(odr.hdr, BCF_DT_ID, vfilt.required_filters[i].c_str()));
    }
  }

  int32_t nVariant = 0, pos = 0;
  int32_t *info_ac = NULL, *tot_ac = NULL;
  int32_t n_ac = 0, n_tac = 0, n_pop = 0;
  char *org_pop = NULL;

  for(int32_t k=0; odr.read(iv); ++k) {

    if (!bcf_is_snp(iv)) {continue;}
    if (iv->pos == pos) {continue;}
    pos = iv->pos;
    // periodic message to user
    if ( k % verbose == 0 )
      notice("Processing %d markers at %s:%d. Recorded %d singletons", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nVariant);

    // unpack FILTER column
    bcf_unpack(iv, BCF_UN_FLT);
    // check --apply-filters
    bool has_filter = req_flt_ids.empty() ? true : false;
    if ( ! has_filter ) {
      for(int32_t i=0; i < iv->d.n_flt; ++i) {
        for(int32_t j=0; j < (int32_t)req_flt_ids.size(); ++j) {
          if ( req_flt_ids[j] == iv->d.flt[i] )
            has_filter = true;
        }
      }
    }
    if ( ! has_filter ) { continue; }
    // check filter logic
    if ( filt != NULL ) {
      int32_t ret = filter_test(filt, iv, NULL);
      if ( filter_logic == FLT_INCLUDE ) { if ( !ret)  has_filter = false; }
      else if ( ret ) { has_filter = false; }
    }
    if ( ! has_filter ) { continue; }

    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    int32_t ac = info_ac[0];
    if (ac > max_ac) {continue;}

    int32_t t_ac = 0;
    if (bcf_get_info_int32(odr.hdr, iv, "TotAC", &tot_ac, &n_tac) >= 0) {
      t_ac = tot_ac[0];
    }

    std::string key;
    if (bcf_get_info_string(odr.hdr, iv, "Pop", &org_pop, &n_pop) < 0) {
      key = "NA";
    }
    else {
      std::string opop(org_pop);
      std::vector<std::string> words;
      std::vector<std::string> pairs;
      std::stringstream ss;
      int32_t x;
      split(words, ",", opop);
      int32_t ct = 0;
      for (auto & v : words) {
        split(pairs, ":", v);
        if (!str2int32(pairs[1], x)) {
          notice("Invalid pop count");
          continue;
        }
        if (x > 0) {
          if (ct > 0) {ss << ",";}
          ss << pairs[0];
          ct++;
        }
      }
      key = ss.str();
    }
    auto ptr = totcount.find(key);
    if (ptr == totcount.end()) {
      totcount[key] = new int32_t[(max_ac+1)*2]{0};
      ptr = totcount.find(key);
    }
    ptr->second[ac*2]++;
    if (t_ac == ac)
      ptr->second[ac*2+1]++;
    nVariant++;
  }

  if (!reg.empty()) {
    out += "_"+reg;
  }
  FILE* wf;

  outf = out + "_tot.count";
  wf = fopen(outf.c_str(), "w");
  for (auto & it : totcount) {
    for ( int32_t i = 1; i < max_ac+1; ++i ) {
      fprintf(wf, "%s\t%d\t%d\t%d\n", it.first.c_str(), i, it.second[i*2], it.second[i*2+1]);
    }
  }
  fclose(wf);

  return 0;
}







