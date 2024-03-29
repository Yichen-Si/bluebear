#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"

#include "pbwt_build.h"

void pbwtSuffix(pbwtCursor& pc, std::vector<bool*>& gtmat, std::vector<int32_t>& posvec, std::vector<int32_t*>& dmat, std::vector<int32_t*>& rmat) {
  int32_t N = (int32_t) posvec.size();
  for (int32_t k = N-1; k >= 0; --k) {
    pc.ForwardsAD_suffix(gtmat[k], posvec[k]);
    memcpy(dmat[k], pc.d, pc.M*sizeof(int32_t));
    pc.ReverseA(rmat[k]);
  }
}

// Input bcf file
// Output suffix pbwt by chunk
int32_t pbwtBuildSuffix(int32_t argc, char** argv) {

  std::string inVcf, reg, out;
  std::string chrom = "chr20";
  int32_t chunksize = 100000;   // Read genotype by chunk
  int32_t store_interval = 5000; // Store snapshot of pbwt
  int32_t nsamples=0, M=0;
  int32_t min_variant = 1;
  int32_t min_ac = 10;
  bool    hyploid = false;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("chr",&chrom, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("chunk-size",&chunksize, "Read, process and write by chunk")
    LONG_INT_PARAM("min-ac",&min_ac, "Minimum allele count to include")
    EXCLUSIVE_PARAM("hyploid",&hyploid, "Assume input is hyploid coded in diploid vcf form and use only the first allele e.g. for chrY")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file prefix")
    LONG_INT_PARAM("store-per-markers",&store_interval, "Snapshot per x markers")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // std::string inVcf = "/net/wonderland/home/ycsi/IBD/ibs0/Data/test/test.10.bcf";
  // std::string reg = "chr20:3100000-3200000";

  int32_t start = 0, end = 300000000;
  std::vector<std::string> v;
  if (!reg.empty()) {
    split(v, ":-", reg);
    if (!str2int32(v[1], start) || !str2int32(v[2], end)) {
      error("Invalid region.");
    }
    chrom = v[0];
  }

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
      break;
    }
  }
  if (nsamples == 0) {
    error("Invalid --region or unreadable --inVcf.");
  }
  M = nsamples;
  if (!hyploid) {
    M *= 2;
  }

  nchunk = ck;
  pbwtCursor pc(M, start + (ck+1) * chunksize + 1);

  notice("Started Reading VCF, identifying %d samples. Will procees in %d chunks and store snapshot every %d markers.", nsamples, nchunk+1, store_interval);

  for (ck = nchunk; ck >= 0; --ck) {

    st = start + ck*chunksize + 1;
    ed = start + ck*chunksize + chunksize;
    if (st == end) continue;
    if (ed > end) ed = end;
    reg = chrom + ":" + std::to_string(st) + "-" + std::to_string(ed);
    notice("Processing %s.", reg.c_str());
    std::vector<GenomeInterval> intervals;
    parse_intervals(intervals, "", reg);
    BCFOrderedReader odr(inVcf, intervals);
    bcf1_t* iv = bcf_init();

    std::vector<bool*> gtmat;
    std::vector<int32_t> positions;

    int32_t *p_gt = NULL;
    int32_t  n_gt = 0;
    int32_t *info_ac = NULL;
    int32_t  n_ac = 0;

    // read marker
    for (int32_t k=0; odr.read(iv); ++k) {
      if (!bcf_is_snp(iv)) {continue;}
      if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
      if (info_ac[0] < min_ac || info_ac[0] > 2*nsamples - min_ac) {continue;}
      // if (k > 10) {break;}
      if (bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
      }

      bool *y;
      y = new bool[M];
      for (int32_t i = 0; i <  nsamples; ++i) {
        int32_t g1 = p_gt[2*i];
        if (hyploid) {
          if (bcf_gt_is_missing(g1)) {
            y[i] = 0;
          } else {
            y[i] = ((bcf_gt_allele(g1) > 0) ? 1 : 0);
          }
        } else {
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
      }
      gtmat.push_back(y);
      positions.push_back(iv->pos+1);
    }

    int32_t N = positions.size();
    if ( N < min_variant ) {
      notice("Observed only %d informative markers. Skipping this chunk...", N);
      continue;
    } else {
      notice("Read %d markers, start extending pbwt for %s.", N, reg.c_str());
    }

    // Build pbwt (suffix) for this block
    // Output by row
    std::string d_outf = out + "_" + reg + "_suffix.dmat";
    std::string a_outf = out + "_" + reg + "_suffix.amat";
    std::ofstream d_wf, a_wf;
    d_wf.open(d_outf);
    a_wf.open(a_outf);
    for (int32_t k = N-1; k >= 0; --k) {
      pc.ForwardsAD_suffix(gtmat[k], positions[k]);

      if (k == N-1 || k == 0 || k % store_interval == 0) {
        d_wf << chrom << '\t' << positions[k];
        for (int32_t j = 0; j < M; ++j)
          d_wf << '\t' << pc.d[j];
        d_wf << '\n';

        a_wf << chrom << '\t' << positions[k];
        for (int32_t j = 0; j < M; ++j)
          a_wf << '\t' << pc.a[j];
        a_wf << '\n';
      }
    }
    d_wf.close();
    a_wf.close();

    for (int32_t k = 0; k < N; ++k) {
      delete [] gtmat[k];
    }
  }


// // Print out test results. Sorted from the first position
// for (uint32_t i = 0; positions[i] <= pc.d[1]; ++i) {
//   int32_t ac = 0;
//   for (int32_t j = 0; j < 10; ++j) {
//     ac += gtmat[i][pc.a[j]];
//   }
//   if (ac == 0) {continue;}
//   std::cout << positions[i]-positions[0] << '\t';
//   for (int32_t j = 0; j < 10; ++j) {
//     std::cout << gtmat[i][pc.a[j]] << ' ';
//   }
//   std::cout << std::endl;
// }
// std::cout << std::endl;
// for (int32_t j = 0; j < M; ++j) {
//   std::cout << pc.d[j]-positions[0] << ' ';
// }
// std::cout << std::endl;

  return 0;
}
