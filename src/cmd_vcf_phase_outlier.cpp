#include "cramore.h"
#include "utils.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "pbwt_build.h"
#include "bp2cm.h"
#include "rare_variant_config.h"

// goal -- For each (common) variant,
//         evaluate how much does it change the haplotype pattern
int32_t HaplotypeChangeAtPt(int32_t argc, char** argv) {
  std::string inVcf, inMap, chrom, reg;
  std::string out;
  int32_t verbose = 10000;
  int32_t min_ac = 10, max_ac = -1, ac_cut = 5;
  int32_t start, end;
  int32_t min_pos = -1, max_pos = -1;
  int32_t chunksize = 1000000;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
    LONG_STRING_PARAM("map",&inMap, "Map file for genetic distance")
    LONG_STRING_PARAM("chr",&chrom, "Chromosome")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("min-pos",&min_pos, "The start point to search for ibd (bp)")
    LONG_INT_PARAM("max-pos",&max_pos, "The end point to search for ibd (bp)")
    LONG_INT_PARAM("ac_cut",&ac_cut, "Ignore AC < ac_cut in haplotype matching")
    LONG_INT_PARAM("min-ac",&min_ac, "Minimum minor allele count to be considered")
    LONG_INT_PARAM("max-ac",&max_ac, "Maximum minor allele count to be considered")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file prefix")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // Genetic map
  bp2cmMap pgmap(inMap, " ", "", -1, -1);
  notice("Read map. min %d; max %d; cst %d; ced %d.", pgmap.minpos, pgmap.maxpos, pgmap.centromere_st, pgmap.centromere_ed);
  if (min_pos < 0)
    min_pos = 0;
  if (max_pos < 0)
    max_pos = pgmap.maxpos;

  // Region to process
  if (reg.empty()) {
    start = min_pos;
    end = max_pos;
  } else {
    std::vector<std::string> v;
    split(v, ":-", reg);
    chrom = v[0];
    if (!str2int32(v[1], start) || !str2int32(v[2], end)) {
      error("Invalid region.");
    }
  }

  // bcf reader
  reg = chrom + ":" + std::to_string(min_pos) + "-" + std::to_string(end+1);
  std::vector<GenomeInterval> intervals;
  parse_intervals(intervals, "", reg);
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();

  // Output pref & suff separately
  std::string out1 = out + "_left.txt";
  std::string out2 = out + "_right.txt";
  FILE* wf;
  wf = fopen(out1.c_str(), "w");
  fprintf(wf, "%s\n", "CHR\tPOS\tAC\tCurrent_Switch_left");

  int32_t nVariant = 0;
  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
  int32_t M = nsamples * 2;

  if (max_ac < 0)
    max_ac = nsamples;

  notice("Started Reading site information from VCF file, identifying %d samples", nsamples);
  if ( nsamples == 0 )
    error("FATAL ERROR: The VCF does not have any samples with genotypes");

  int32_t* p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t *info_ac = NULL;
  int32_t  n_ac = 0;
  int32_t ac = 0;

  // Build pbwt from beginning.
  pbwtCursor pc(M, min_pos);

  // forward pbwt
  for (int32_t k=0; odr.read(iv); ++k) {
    if ( k % verbose == 0 )
      notice("Forward pbwt. Processing %d markers at %s:%d. Recorded %d variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nVariant);
    if (bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    ac = info_ac[0];
    if (ac > nsamples)
      ac = M - ac;

    bool y[M];
    for (int32_t i = 0; i <  nsamples; ++i) {
      int32_t g1 = p_gt[2*i];
      int32_t g2 = p_gt[2*i+1];
      if (bcf_gt_is_missing(g1)) {
        y[2*i] = 0;
      } else {
        y[2*i] = ((bcf_gt_allele(g1) > 0) ? 1 : 0);
      }
      if (bcf_gt_is_missing(g2)) {
        y[2*i+1] = 0;
      } else {
        y[2*i+1] = ((bcf_gt_allele(g2) > 0) ? 1 : 0);
      }
    }

    if ( ac >= min_ac && ac <= max_ac ) {
      // Basedd on haplotype matching order before this position
      int32_t pre_allele = y[pc.a[0]];
      int32_t switch_ct = 0;
      for (int32_t i = 1; i < M; ++i) {
        if (y[pc.a[i]] != pre_allele) {
          switch_ct += 1;
          pre_allele = y[pc.a[i]];
        }
      }
      fprintf(wf,"%s\t%d\t%d\t%d\n", bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, ac, switch_ct);
      nVariant++;
    }
    // Update pbwt cursor
    if (ac < ac_cut) {continue;}
    pc.ForwardsAD_prefix(y, iv->pos+1);
  }
  fclose(wf);

  wf = fopen(out2.c_str(), "w");
  fprintf(wf, "%s\n", "CHR\tPOS\tAC\tCurrent_Switch_right");

  // backward pbwt
  pc.Reset(M, max_pos);

  int32_t ed = max_pos;
  int32_t st = ed - chunksize + 1;

  // Have to read by chunk
  while ( ed > start ) {
    if (st < start) st = start;
    reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
    notice("Processing %s.", reg.c_str());
    parse_intervals(intervals, "", reg);
    odr.jump_to_interval(intervals[0]);
    bcf_clear(iv);

    std::vector<bool*> gtmat;
    std::vector<int32_t> positions;
    std::vector<int32_t> mac;

    int32_t *p_gt = NULL;
    int32_t  n_gt = 0;
    int32_t *info_ac = NULL;
    int32_t n_ac = 0;

    // read marker
    for (int32_t k=0; odr.read(iv); ++k) {
      if (bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
      }
      if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
      ac = info_ac[0];
      if (ac > nsamples)
        ac = M - ac;
      if (ac < ac_cut) {continue;}
      bool *y;
      y = new bool[M];
      for (int32_t i = 0; i <  nsamples; ++i) {
        int32_t g1 = p_gt[2*i];
        int32_t g2 = p_gt[2*i+1];
        if (bcf_gt_is_missing(g1)) {
          y[2*i] = 0;
        } else {
          y[2*i] = ((bcf_gt_allele(g1) > 0) ? 1 : 0);
        }
        if (bcf_gt_is_missing(g2)) {
          y[2*i+1] = 0;
        } else {
          y[2*i+1] = ((bcf_gt_allele(g2) > 0) ? 1 : 0);
        }
      }
      gtmat.push_back(y);
      positions.push_back(iv->pos+1);
      mac.push_back(ac);
    }

    int32_t N = positions.size();
    if ( N < 1 ) {
      notice("Observed 0 informative markers. Skipping this chunk...");
      ed -=chunksize;
      st = ed  - chunksize + 1;
      continue;
    } else {
      notice("Read %d markers, start backward pbwt for %s.", N, reg.c_str());
    }

    for (int32_t k = N-1; k >= 0; --k) {
      if ( positions[k] < end && mac[k] >= min_ac && mac[k] <= max_ac ) {
        // Basedd on haplotype matching order before this position
        int32_t pre_allele = gtmat[k][pc.a[0]];
        int32_t switch_ct = 0;
        for (int32_t i = 1; i < M; ++i) {
          if (gtmat[k][pc.a[i]] != pre_allele) {
            switch_ct += 1;
            pre_allele = gtmat[k][pc.a[i]];
          }
        }
        fprintf(wf,"%s\t%d\t%d\t%d\t\n", bcf_hdr_id2name(odr.hdr, iv->rid), positions[k], mac[k], switch_ct);
        nVariant++;
      }
      pc.ForwardsAD_suffix(gtmat[k], positions[k]);
    } // Finish building pbwt for this block
    for (int32_t k = 0; k < N; ++k) {
      delete [] gtmat[k];
    }
    ed -=chunksize;
    st = ed  - chunksize + 1;
  }
  fclose(wf);

  return 0;
}
