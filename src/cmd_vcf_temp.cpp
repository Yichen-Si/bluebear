#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
// #include "bcf_ordered_writer.h"
#include "hts_utils.h"
#include <bits/stdc++.h>
// #include "ibs0.h"
#include <sstream>


int32_t temp(int32_t argc, char** argv) {

  std::string inVcf, reg, out;
  int32_t max_ac=2, min_ac=2;
  int32_t verbose = 10000;

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("max-ac",&max_ac, "Maximum sample allele count to consider")
    LONG_INT_PARAM("min-ac",&min_ac, "Minimum sample allele count to consider")
    LONG_INT_PARAM("verbose",&verbose, "Periodic report when scanning VCF")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if (inVcf.empty() || out.empty()) {
    error("[E:%s:%d %s] --in-vcf, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // output
  FILE * wf;
  wf = fopen(out.c_str(), "w");
  fprintf(wf, "%s\n", "CHR\tPOS\tAC\tID1\tID2\tRefAC");

  // bcf reader
  std::vector<GenomeInterval> intervals;
  parse_intervals(intervals, "", reg);
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();
  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
  char** id_ptr = bcf_hdr_get_samples(odr.hdr);
  std::vector<std::string> id_samples(id_ptr, id_ptr+nsamples);

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

  int32_t* p_gt = NULL;
  int32_t* info_ac = NULL;
  int32_t* info_an = NULL;
  int32_t* ref_ac = NULL;
  int32_t n_gt = 0, n_ac = 0, n_an = 0, n_ref = 0, pos = 0;
  int32_t nVariant = 0;

  for(int32_t k=0; odr.read(iv); ++k) {
    if (iv->pos+1 == pos) {
      continue; // Jump over tri-allelic sites
    }
    // periodic message to user
    if ( k % verbose == 0 )
      notice("Processing %d markers at %s:%d. Recorded %d rare variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nVariant);

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
    // extract genotype and apply genotype level filter
    if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (info_ac[0] < 2) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
    int32_t ac = 0, an = 0;
    ac = info_ac[0]; an = info_an[0];
    if (ac > an/2) {
      ac = an - ac;
    }
    if (ac < min_ac || ac > max_ac) {continue;}

    std::vector<std::string> carry;
    for(int32_t i=0; i < nsamples; ++i) {
      int32_t g1 = p_gt[2*i];
      int32_t g2 = p_gt[2*i+1];
      int32_t geno;
      if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
        continue;
      } else {
        geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
        if (geno == 1) {
          carry.push_back( id_samples[i] );
        }
      }
    } // Found carriers (have to be het)
    // TODO: here we ignore rare varaiants with any hom carrier
    if ((int32_t) carry.size() != ac) {continue;}

    pos = iv->pos+1;
    int32_t ref = -1;
    if (bcf_get_info_int32(odr.hdr, iv, "full_data_AC", &ref_ac, &n_ref) >= 0) {
      ref = ref_ac[0];
    }

      // Find pairwise ibs0 w/in a short limit
    for (int32_t i = 0; i < ac-1; ++i) {
      for (int32_t j = i+1; j < ac; ++j) {
        fprintf(wf, "%s\t%d\t%d\t%s\t%s\t%d\n", bcf_hdr_id2name(odr.hdr, iv->rid), pos, ac, carry[i].c_str(), carry[j].c_str(), ref);
      }
    }
    nVariant++;
  }

  fclose(wf);
  return 0;
}














