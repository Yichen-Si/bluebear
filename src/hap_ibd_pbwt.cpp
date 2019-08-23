#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"

#include "pbwt_build.h"

// Input bcf file
// Matched haplotype length preceeding a shared rare variant
int32_t hapIBDpbwtLeft(int32_t argc, char** argv) {

  std::string inVcf, reg, out, outf;
  int32_t verbose = 10000;
  int32_t max_ac  = 3;     // Max AC for anchor IBD.
  int32_t min_ac  = 4;     // Ignore the phase of AC <= min_AC.
  int32_t nsamples=0, M=0;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("max-ac",&max_ac, "Max. AC as anchor for IBD")
    LONG_INT_PARAM("min-ac",&min_ac, "Min. AC to consider in haplotype matching")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file prefix")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // bcf reader
  std::vector<GenomeInterval> intervals;
  if ( !reg.empty() ) {
    parse_intervals(intervals, "", reg);
  }
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();
  nsamples = bcf_hdr_nsamples(odr.hdr);
  M = nsamples * 2;

  // Initialize pbwt Cursor
  pbwtCursor pc(M, 1);

  // Set up output
  outf = out + ".left.list";
  htsFile* wf = hts_open(outf.c_str(), "w");

  int32_t *p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t *info_ac = NULL;
  int32_t *info_an = NULL;
  int32_t n_ac = 0, n_an = 0;
  bool *y;
  y = new bool[M];

  notice("Started Reading VCF, identifying %d samples.", nsamples);

  // read marker
  for (int32_t k=0; odr.read(iv); ++k) {

    if (bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
    n_ac = info_ac[0];
    n_an = info_an[0];
    bool ifflip = 0;
    if (n_ac > n_an/2) {
      ifflip = 1;
      n_ac = n_an - n_ac;
    }

    if (n_ac <= max_ac) { // Check haplotype matching before this position
      std::vector<int32_t> idlist;
      std::vector<int32_t> hapcarry;
      for (int32_t i = 0; i <  nsamples; ++i) {
        int32_t g1 = p_gt[2*i];
        int32_t g2 = p_gt[2*i+1];
        if (!bcf_gt_is_missing(g1) && !bcf_gt_is_missing(g2)) {
          if (ifflip) {
            if (bcf_gt_allele(g1) == 0 || bcf_gt_allele(g2) == 0) {
              idlist.push_back(i);
              hapcarry.push_back((bcf_gt_allele(g1) == 0)? 0 : 1);
            }
          } else {
            if (bcf_gt_allele(g1) > 0 || bcf_gt_allele(g2) > 0) {
              idlist.push_back(i);
              hapcarry.push_back((bcf_gt_allele(g1) > 0)? 0 : 1);
            }
          }
        }
      }
      if ((int32_t) idlist.size() < n_ac) {continue;}
      int32_t rvec[M];
      pc.ReverseA(rvec);
      int32_t h11,h12,h21,h22;
      int32_t d11,d12,d21,d22;
      for (int32_t i = 0; i < n_ac-1; i++) {
        h11 = rvec[i * 2 + hapcarry[i]];
        h12 = rvec[i * 2 + 1 - hapcarry[i]];
        for (int32_t j = i+1; j < n_ac; j++) {
          h21 = rvec[j * 2 + hapcarry[j]];
          h22 = rvec[j * 2 + 1 - hapcarry[j]];
          d11 = iv->pos+1 - pc.Dist_pref(h11,h21);
          d12 = iv->pos+1 - pc.Dist_pref(h11,h22);
          d21 = iv->pos+1 - pc.Dist_pref(h12,h21);
          d22 = iv->pos+1 - pc.Dist_pref(h12,h22);
          hprintf(wf,"%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\n", iv->pos+1, n_ac, odr.hdr->id[BCF_DT_SAMPLE][idlist[i]].key, odr.hdr->id[BCF_DT_SAMPLE][idlist[j]].key, d11,d12,d21,d22);
        }
      }
    }
    if (n_ac < min_ac) {continue;}

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
    if (ifflip) {
      for (int32_t i = 0; i < 2*nsamples; ++i) {
        y[i] = !y[i];
      }
    }
    pc.ForwardsAD_prefix(y, iv->pos+1);
    if ( k % verbose == 0 )
    notice("Processing %d markers at %s:%d.", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
  }

  delete [] y;
  hts_close(wf);

  return 0;
}




// Input bcf file
// Matched haplotype length succeeding a shared rare variant
int32_t hapIBDpbwtRight(int32_t argc, char** argv) {

  std::string inVcf, reg, out, outf, chrom;
  int32_t verbose = 10000;
  int32_t max_ac  = 3;     // Max AC for anchor IBD.
  int32_t min_ac  = 3;     // Ignore the phase of AC <= min_AC.
  int32_t nsamples=0, M=0;
  int32_t chunksize = 500000;
  int32_t start = 0, end = 65000000;
  int32_t min_variant = 1;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("max-ac",&max_ac, "Max. AC as anchor for IBD")
    LONG_INT_PARAM("min-ac",&min_ac, "Min. AC to consider in haplotype matching")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file prefix")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() || reg.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --out, --reg are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // Initial info
  std::vector<std::string> v;
  split(v, ":-", reg);
  if (!str2int32(v[1], start) || !str2int32(v[2], end)) {
    error("Invalid region.");
  }
  chrom = v[0];

  int32_t nchunk = (end-start)/chunksize;
  int32_t ck=nchunk;
  int32_t st, ed;

  for (ck = nchunk; ck >= 0; --ck) {
    st = start + ck*chunksize + 1;
    ed = start + ck*chunksize + chunksize;
    reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
    std::vector<GenomeInterval> intervals;
    parse_intervals(intervals, "", reg);
    BCFOrderedReader odr(inVcf, intervals);
    bcf1_t* iv = bcf_init();
    if (odr.read(iv)) {
      nsamples = bcf_hdr_nsamples(odr.hdr);
      M = nsamples * 2;
      break;
    }
  }
  if (M == 0) {
    error("No individual is read from --inVcf.");
  }
  nchunk = ck + 1;

  // Initialize pbwt Cursor
  pbwtCursor pc(M, end);

  // Set up output
  outf = out + ".right.list";
  htsFile* wf = hts_open(outf.c_str(), "w");

  // Have to read by chunk
  for (ck = 0; ck <= nchunk; ++ck) {

    ed = end - ck*chunksize;
    st = ed  - chunksize + 1;
    if (ed <= start) {break;}
    if (st < start) st = start;
    reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
    notice("Processing %s.", reg.c_str());
    std::vector<GenomeInterval> intervals;
    parse_intervals(intervals, "", reg);
    BCFOrderedReader odr(inVcf, intervals);
    bcf1_t* iv = bcf_init();

    std::vector<bool*> gtmat;
    std::vector<int32_t> positions;
    std::vector<int32_t> mac;

    int32_t *p_gt = NULL;
    int32_t  n_gt = 0;
    int32_t *info_ac = NULL;
    int32_t *info_an = NULL;
    int32_t n_ac = 0, n_an = 0;

    // read marker
    for (int32_t k=0; odr.read(iv); ++k) {
      if (bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
      }
      if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
      if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
      n_ac = info_ac[0];
      n_an = info_an[0];
      bool ifflip = 0;
      if (n_ac > n_an/2) {
        ifflip = 1;
        n_ac = n_an - n_ac;
      }
      if (n_ac < 2) {continue;}

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
      if (ifflip) {
        for (int32_t i = 0; i <  2*nsamples; ++i) {
          y[i] = !y[i];
        }
      }

      gtmat.push_back(y);
      positions.push_back(iv->pos+1);
      mac.push_back(n_ac);
    }
    int32_t N = positions.size();
    if ( N < min_variant ) {
      notice("Observed only %d informative markers. Skipping this chunk...", N);
      continue;
    } else {
      notice("Read %d markers, start extending pbwt for %s.", N, reg.c_str());
    }

    for (int32_t k = N-1; k >= 0; --k) {

      if (mac[k] <= max_ac) { // Check haplotype matching starting from this position
        std::vector<int32_t> idlist;
        std::vector<int32_t> hapcarry;
        for (int32_t i = 0; i <  nsamples; ++i) {
          if (gtmat[k][2*i] + gtmat[k][2*i+1] == 1) {
            idlist.push_back(i);
            hapcarry.push_back((gtmat[k][2*i] > 0) ? 0 : 1);
          }
        }
        if ((int32_t) idlist.size() < mac[k]) {continue;}
        int32_t rvec[M];
        pc.ReverseA(rvec);
        int32_t h11,h12,h21,h22;
        int32_t d11,d12,d21,d22;
        for (int32_t i = 0; i < mac[k]-1; i++) {
          h11 = rvec[i * 2 + hapcarry[i]];
          h12 = rvec[h11 + 1 - 2 * (h11 % 2)];
          for (int32_t j = i+1; j < mac[k]; j++) {
            h21 = rvec[j * 2 + hapcarry[j]];
            h22 = rvec[h21 + 1 - 2 * (h21 % 2)];
            d11 = pc.Dist_suff(h11,h21) - positions[k];
            d12 = pc.Dist_suff(h11,h22) - positions[k];
            d21 = pc.Dist_suff(h12,h21) - positions[k];
            d22 = pc.Dist_suff(h12,h22) - positions[k];
            hprintf(wf,"%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\n", positions[k], mac[k], odr.hdr->id[BCF_DT_SAMPLE][idlist[i]].key, odr.hdr->id[BCF_DT_SAMPLE][idlist[j]].key, d11,d12,d21,d22);
          }
        }
      }
      if (mac[k] < min_ac) {continue;}

      pc.ForwardsAD_suffix(gtmat[k], positions[k]);

    } // Finish building pbwt for this block

    for (int32_t k = 0; k < N; ++k) {
      delete [] gtmat[k];
    }

  }

  hts_close(wf);
  return 0;
}















