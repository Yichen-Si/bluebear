#include "cramore.h"
#include "utils.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "ibs0.h"
#include "rare_variant_ibs0.h"
#include <stdio.h>
#include <iomanip>
#include <sstream>

// goal -- Get no-ibs0 length between carriers of different alleles
//         at a tri-allelic site
int32_t IBS0Triallelic(int32_t argc, char** argv) {
  std::string inVcf, inMap, chrom, reg;
  std::string out;
  int32_t min_hom_gts = 1;
  int32_t verbose = 10000;
  int32_t start = -1, end = -1;
  int32_t window_size = 500000;
  int32_t cst = -1, ced = -1;
  int32_t ck_len = 500000;
  int32_t bp_limit = 2000000;
  int32_t max_pair = 3;

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input data", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("chrom",&chrom, "Chromosome to sample ibs0 length")
    LONG_STRING_PARAM("map",&inMap, "Genetic map for genetic distance")
    LONG_STRING_PARAM("reg",&reg, "Region to work on")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")
    LONG_INT_PARAM("centromere-st",&cst, "Start position of centromere")
    LONG_INT_PARAM("centromere-ed",&ced, "End position of centromere")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("min-hom",&min_hom_gts, "Minimum number of homozygous genotypes to be counted for IBS0")
    LONG_INT_PARAM("bp-limit",&bp_limit, "Max distance to look for IBS0 (bp)")
    LONG_INT_PARAM("window-size",&window_size, "Length of the sliding window to procee (bp)")
    LONG_INT_PARAM("max-pair",&max_pair, "Maxiimum number of pairs to record for each multi-allelic site")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || inMap.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --inMap, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }


  // Genetic map
  bp2cmMap pgmap(inMap, " ", "", cst, ced);
  notice("Read map file and detect centromere: %d-%d", pgmap.centromere_st, pgmap.centromere_ed);

  // Region to process
  if (!reg.empty()) {
    std::vector<std::string> v;
    split(v, ":-", reg);
    chrom = v[0];
    if (!str2int32(v[1], start) || !str2int32(v[2], end)) {
      error("Invalid region.");
    }
  }

  // Input vcf
  std::vector<GenomeInterval> intervals;
  parse_intervals(intervals, "", reg);
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();
  bcf1_t* iv_pre = bcf_init();

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

  // Output
  FILE * wf;
  wf = fopen(out.c_str(), "w");
  fputs("CHR\tPOS\tAlleles\tAC1\tAC2\tID1\tID2\tIBS0\tOverlap\n", wf);

  // May need to randomly sample individual pairs
  std::random_device rd;
  std::mt19937 rng(rd());
  std::set<int32_t> chosen;
  std::vector<int32_t> pair(2);
  int32_t leftibs0, rightibs0;


  if (start < 0)
    start = pgmap.minpos;
  end = start + window_size;
  reg = chrom + ":" + std::to_string(start) + "-" + std::to_string(end);

  // IBS0 reader
  IBS0lookup ibs0finder(inVcf, reg, pgmap, bp_limit, ck_len, min_hom_gts);
  if (ibs0finder.start_que.size() < 1) {
    notice("Not enough variants in the given region. Stopped without output.");
    return 0;
  }
  notice("Set up ibs0finder for region %s", reg.c_str());

  // Read first SNP
  while(odr.read(iv_pre)) {
     // unpack FILTER column
    bcf_unpack(iv_pre, BCF_UN_FLT);
    // check --apply-filters
    bool has_filter = req_flt_ids.empty() ? true : false;
    if ( ! has_filter ) {
      for(int32_t i=0; i < iv_pre->d.n_flt; ++i) {
        for(int32_t j=0; j < (int32_t)req_flt_ids.size(); ++j) {
          if ( req_flt_ids[j] == iv_pre->d.flt[i] )
            has_filter = true;
        }
      }
    }
    if ( ! has_filter ) { continue; }
    // check filter logic
    if ( filt != NULL ) {
      int32_t ret = filter_test(filt, iv_pre, NULL);
      if ( filter_logic == FLT_INCLUDE ) { if ( !ret)  has_filter = false; }
      else if ( ret ) { has_filter = false; }
    }
    if ( ! has_filter ) { continue; }
    if (!bcf_is_snp(iv_pre)) { continue; }
    break;
  }
  int32_t pos = iv_pre->pos;

  for(int32_t k=0; odr.read(iv); ++k) {
    // periodic message to user
    if ( k % verbose == 0 )
      notice("Processing %d markers at %s:%d. Recorded %d tri-allelic site", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nVariant);

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
    if (!bcf_is_snp(iv)) { continue; }

    if (iv->pos == pos) {
      // Found one tri-allelic site

      int32_t* p_gt = NULL;
      int32_t  n_gt = 0;
      int32_t *info_ac = NULL;
      int32_t *info_an = NULL;
      int32_t n_ac = 0, n_an = 0;

      int32_t minor1 = 1, minor2 = 1;
      int32_t ac1 = 0, ac2 = 0, an = 0;

      // Extract carriers of allele 1
      if ( bcf_get_genotypes(odr.hdr, iv_pre, &p_gt, &n_gt) < 0 ) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv_pre->rid), iv_pre->pos+1);
      }
      if (bcf_get_info_int32(odr.hdr, iv_pre, "AC", &info_ac, &n_ac) < 0) {continue;}
      if (bcf_get_info_int32(odr.hdr, iv_pre, "AN", &info_an, &n_an) < 0) {continue;}
      ac1 = info_ac[0]; an = info_an[0];
      if (ac1 > an/2) {
        ac1 = an - ac1;
        minor1 = 0;
      }
      std::vector<int32_t> carry1;
      for(int32_t i=0; i < nsamples; ++i) {
        int32_t g1 = p_gt[2*i];
        int32_t g2 = p_gt[2*i+1];
        int32_t geno;
        if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
          continue;
        } else {
          geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
          if (geno == minor1 || geno == minor1*2) {
            carry1.push_back(i);
          }
        }
      } // Found carriers from allele 1
      std::set<int32_t> cset1(carry1.begin(), carry1.end());

      p_gt = NULL; info_ac = NULL; info_an = NULL;
      n_gt = 0; n_ac = 0; n_an = 0;

      // Extract genotype from allele 2
      if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
      }
      if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
      if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
      ac2 = info_ac[0]; an = info_an[0];
      if (ac2 > an/2) {
        ac2 = an - ac2;
        minor2 = 0;
      }
      std::vector<int32_t> carry2, ovl;
      for(int32_t i=0; i < nsamples; ++i) {
        int32_t g1 = p_gt[2*i];
        int32_t g2 = p_gt[2*i+1];
        int32_t geno;
        if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
          continue;
        } else {
          geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
          if (geno == minor2 || geno == minor2*2) {
            if (cset1.find(i) == cset1.end()) {
              carry2.push_back(i);
            } else {
              ovl.push_back(i);
            }
          }
        }
      } // Found carriers from allele 2

      // Decide which pairs to get IBS0
      if (carry1.size()>0 && carry2.size()>0) {
        std::set<int32_t> pair1, pair2;
        int32_t n_pair = max_pair;
        n_pair = std::min(n_pair, (int32_t) carry1.size());
        n_pair = std::min(n_pair, (int32_t) carry2.size());

        NchooseK(carry1.size(), n_pair, pair1, rng, 0);
        NchooseK(carry2.size(), n_pair, pair2, rng, 0);

        std::vector<int32_t> pair1v(n_pair), pair2v(n_pair);
        int32_t it = 0;
        for (auto v : pair1) {
          pair1v[it]=carry1[v]; it++;
        }
        it = 0;
        for (auto v : pair2) {
          pair2v[it]=carry2[v]; it++;
        }

        std::vector<int32_t> bplist(n_pair);
        std::vector<double> cmlist(n_pair);
        for (it=0; it<n_pair; ++it) {
          leftibs0 = ibs0finder.FindIBS0(pair1v[it],pair2v[it],iv->pos+1,1);
          rightibs0= ibs0finder.FindIBS0(pair1v[it],pair2v[it],iv->pos+1,0);
          if (leftibs0 <= 0 || rightibs0 <= 0) {
            bplist[it] = bp_limit*2+window_size;
            cmlist[it] = pgmap.bpinterval2cm(pos-bp_limit,iv->pos+bp_limit);
          } else {
            bplist[it] = rightibs0 - leftibs0;
            cmlist[it] = pgmap.bpinterval2cm(leftibs0,rightibs0);
          }
        }

        std::string ref,alt1,alt2;
        if (minor2) {
          ref = iv->d.allele[0];
          alt2= iv->d.allele[1];
        } else {
          ref = iv->d.allele[1];
          alt2= iv->d.allele[0];
        }
        if (minor1) {
          alt1= iv_pre->d.allele[1];
        } else {
          alt1= iv_pre->d.allele[0];
        }
        std::string allele = ref + "," + alt1 + "," + alt2;

        std::stringstream ss;
        for (it=0; it<n_pair; ++it) {
          ss << pair1v[it] << "," << pair2v[it] << ":" << bplist[it] << "," << std::fixed << std::setprecision(4) << cmlist[it] << ";";
        }

        std::stringstream id1, id2;
        id1 << carry1[0];
        it = 1;
        while(it < (int32_t) carry1.size()) {
          id1 << "," << carry1[it];
          it++;
        }
        id2 << carry2[0];
        it = 1;
        while(it < (int32_t) carry2.size()) {
          id2 << "," << carry2[it];
          it++;
        }

        std::stringstream overlap;
        if (ovl.size() == 0) {overlap << 0;}
        else {
          for (auto v : ovl) {
            overlap << v << ";";
          }
        }

// std::cout << iv->pos+1 << '\t' << allele << '\t' << ac1 << '\t' << ac2 << '\t' << ss.str() << '\t' << overlap.str() << '\n';

        fprintf(wf, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n", chrom.c_str(), iv->pos+1, allele.c_str(), ac1, ac2, id1.str().c_str(), id2.str().c_str(), ss.str().c_str(), overlap.str().c_str());
        nVariant++;
      }

      pos = -1;
    } else {
      pos = iv->pos;
      bcf_destroy(iv_pre);
      iv_pre = bcf_dup(iv);
    }

    if (iv->pos > end) { // Update
      start = end + 1;
      end = start + window_size;
      while (pgmap.overlap_centro(start,end) > 0 && end < pgmap.maxpos) {
        start = end + 1;
        end = start + window_size;
      }
      if (start > pgmap.maxpos) {break;}
      reg = chrom + ":" + std::to_string(start) + "-" + std::to_string(end);
      ibs0finder.Update(reg);
      notice("Update ibs0finder for region %s", reg.c_str());
    }

  }

  pgmap.clear();
  fclose(wf);

  return 0;
}

