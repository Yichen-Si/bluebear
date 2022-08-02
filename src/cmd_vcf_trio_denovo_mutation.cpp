#include <stdio.h>
#include <iostream>
#include <fstream>

#include "cramore.h"
#include "utils.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"

struct trio {
    std::string family_id;
    int32_t paternal_id, maternal_id, child_id;
    int32_t denovo_total;
};

int32_t trioDenovo(int32_t argc, char** argv) {
  std::string inVcf, inPed;
  std::string out, outf;
  int32_t verbose = 100000;
  std::string reg;
  int32_t n_variant = 0, n_family = 0;

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file of genotype")
    LONG_STRING_PARAM("in-ped",&inPed, "Input VCF/BCF file of family info")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file prefix")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() || inPed.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --in-ped and --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // bcf reader
  std::vector<GenomeInterval> intervals;
  if ( !reg.empty() ) {
    parse_intervals(intervals, "", reg);
  }
  BCFOrderedReader odr_c(inVcf, intervals);
  bcf1_t* iv_c = bcf_init();

  int32_t n_samples = bcf_hdr_nsamples(odr_c.hdr);
  char** id_ptr = bcf_hdr_get_samples(odr_c.hdr);
  std::vector<std::string> id(id_ptr, id_ptr+n_samples);
  std::map<std::string, int32_t> id_map;
  for (int32_t i = 0; i < n_samples; ++i) {
      id_map[id[i]] = i;
  }
  if ( n_samples == 0 )
    error("FATAL ERROR: The VCF does not have any samples with genotypes");

  // family info
  std::vector<trio> trios;
  std::ifstream rf;
  std::string line;
  std::vector<std::string> v;
  rf.open(inPed, std::ifstream::in);
  while (std::getline(rf, line)) {
      split(v, "\t", line);
      if (v.size() < 4) {
          error("Invalid input family information.");
      }
      auto p_it = id_map.find(v[2]);
      if (p_it == id_map.end()) {
          continue;
      }
      auto c_it = id_map.find(v[1]);
      auto m_it = id_map.find(v[3]);
      if (c_it == id_map.end() || m_it == id_map.end()) {
          continue;
      }
      trio rec = {v[0], p_it->second, m_it->second, c_it->second, 0};
      trios.push_back(rec);
      n_family++;
  }
  rf.close();
  if ( n_family == 0 )
    error("FATAL ERROR: Input does not contain any trio");

  // Output
  outf = out + ".bed";
  FILE *wf;
  wf = fopen(outf.c_str(), "w");
  fprintf(wf, "#chrom\tchromStart\tchromEnd\tREF\tALT\tAC\tFamily_id\tDenovo\tPaternal_gt\tMaternal_gt\tChild_gt\n");

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
  if ( filter_logic != 0 ) {
    filter_init(odr_c.hdr, filter_str.c_str());
  }

  // handle --apply-filtrs
  std::vector<int32_t> req_flt_ids;
  if ( !vfilt.required_filters.empty() ) {
    for(int32_t i=0; i < (int32_t)vfilt.required_filters.size(); ++i) {
      req_flt_ids.push_back(bcf_hdr_id2int(odr_c.hdr, BCF_DT_ID, vfilt.required_filters[i].c_str()));
    }
  }

  notice("Started Reading site information from VCF file, identifying %d samples and %d trios", n_samples, n_family);

  int32_t *p_gt = NULL, *p_ac = NULL, *p_an = NULL;
  int32_t  n_gt = 0, n_ac = 0, n_an = 0, ac = 0, an = 0;
  int32_t  n_denovo = 0;

  for(int32_t k=0; odr_c.read(iv_c); ++k) {

    if ( k % verbose == 0 )
      notice("Processing vcf: %d markers at %s:%d. Found %d de novo incidences", n_variant, bcf_hdr_id2name(odr_c.hdr, iv_c->rid), iv_c->pos+1, n_denovo);
    bcf_unpack(iv_c, BCF_UN_FLT);
    if (!bcf_is_snp(iv_c)) { continue; }
    // check --apply-filters
    bool has_filter = req_flt_ids.empty() ? true : false;
    if ( ! has_filter ) {
      //notice("%d %d", iv_c->d.n_flt, (int32_t)req_flt_ids.size());
      for(int32_t i=0; i < iv_c->d.n_flt; ++i) {
        for(int32_t j=0; j < (int32_t)req_flt_ids.size(); ++j) {
          if ( req_flt_ids[j] == iv_c->d.flt[i] )
            has_filter = true;
        }
      }
    }
    if ( ! has_filter ) { continue; }
    // check filter logic
    if ( filt != NULL ) {
      int32_t ret = filter_test(filt, iv_c, NULL);
      if ( filter_logic == FLT_INCLUDE ) { if ( !ret)  has_filter = false; }
      else if ( ret ) { has_filter = false; }
    }
    if ( ! has_filter ) { continue; }

    if (bcf_get_info_int32(odr_c.hdr, iv_c, "AC", &p_ac, &n_ac) < 0) {continue;}
    if (bcf_get_info_int32(odr_c.hdr, iv_c, "AN", &p_an, &n_an) < 0) {continue;}
    ac = p_ac[0];
    an = p_an[0];
    if (ac == 0) { continue; }
    // if (an < n_samples) { continue; }
    // extract genotype and apply genotype level filter
    if ( bcf_get_genotypes(odr_c.hdr, iv_c, &p_gt, &n_gt) < 0 ) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr_c.hdr, iv_c->rid), iv_c->pos+1);
    }
    for (int32_t i = 0; i < n_family; ++i) {
        std::vector<int32_t> gtv(6, 0);
        std::vector<int32_t> dgt(3, 0);
        int32_t alt_tot = 0;
        trio& rec = trios[i];
        gtv[0]=p_gt[2*rec.child_id];
        gtv[1]=p_gt[2*rec.child_id+1];
        gtv[2]=p_gt[2*rec.paternal_id];
        gtv[3]=p_gt[2*rec.paternal_id+1];
        gtv[4]=p_gt[2*rec.maternal_id];
        gtv[5]=p_gt[2*rec.maternal_id+1];
        for (int32_t j = 0; j < 3; ++j) {
            if ( bcf_gt_is_missing(gtv[j*2]) || bcf_gt_is_missing(gtv[j*2+1])) {
                // If any allele is missing, don't make decision
            }
            else {
                gtv[j*2] = (bcf_gt_allele(gtv[j*2]) > 0) ? 1 : 0;
                gtv[j*2+1] = (bcf_gt_allele(gtv[j*2+1]) > 0) ? 1 : 0;
                dgt[j] = gtv[j*2] + gtv[j*2+1];
                alt_tot += dgt[j];
            }
        }
        if (alt_tot == 0) {continue;}
        int32_t delta = 0;
        if (dgt[1] == 0 && dgt[2] == 0 && dgt[0] > 0) {
            delta = dgt[0];
        }
        if (dgt[1] + dgt[2] == 1 && dgt[0] == 2) {
            delta = 1;
        }
        if (dgt[1] + dgt[2] == 3 && dgt[0] == 0) {
            delta = -1;
        }
        if (dgt[1] + dgt[2] == 4 && dgt[0] < 2) {
            delta = dgt[0]-2;
        }
        if (delta != 0) {
            char** alleles = bcf_get_allele(iv_c);
            fprintf(wf, "%s\t%d\t%d\t%s\t%s\t%d\t%s\t%d\t%d\t%d\t%d\n", bcf_hdr_id2name(odr_c.hdr, iv_c->rid), iv_c->pos, iv_c->pos+1, alleles[0], alleles[1], ac, rec.family_id.c_str(), delta, dgt[1], dgt[2], dgt[0]);
            rec.denovo_total++;
            n_denovo++;
        }
    }
    n_variant++;
  }
  if (p_gt) free(p_gt);
  if (n_ac) free(p_ac);
  // if (n_an) free(p_an);
  fclose (wf);

  // Output family summary
  outf = out + ".trio_summary.tsv";
  wf = fopen(outf.c_str(), "w");
  fprintf(wf, "Family_id\tChild_id\tPaternal_id\tMaternal_id\tDenovo_site\n");
  for (int32_t i = 0; i < n_family; ++i) {
      trio& rec = trios[i];
      fprintf(wf, "%s\t%s\t%s\t%s\t%d\n", rec.family_id.c_str(),
      id[rec.child_id].c_str(), id[rec.paternal_id].c_str(), id[rec.maternal_id].c_str(),
      rec.denovo_total);
  }
  fclose (wf);
  notice("Finish processing %d markers. Found %d de novo incidences", n_variant, n_denovo);
  return 0;
}
