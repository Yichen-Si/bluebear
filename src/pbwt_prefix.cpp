#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"

#include "pbwt_build.h"

void pbwtPrefix(pbwtCursor& pc, std::vector<bool*>& gtmat, std::vector<int32_t>& posvec, std::vector<int32_t*>& dmat, std::vector<int32_t*>& rmat) {
  int32_t N = (int32_t) posvec.size();
  for (int32_t k = 0; k < N; ++k) {
    pc.ForwardsAD_prefix(gtmat[k], posvec[k]);
    memcpy(dmat[k], pc.d, pc.M*sizeof(int32_t));
    pc.ReverseA(rmat[k]);
  }
}

// Input bcf file
// Output snapshot of prefix pbwt
int32_t pbwtBuildPrefix(int32_t argc, char** argv) {

  std::string inVcf, reg, out;
  std::string chrom = "chr20";
  int32_t verbose = 10000;
  int32_t store_interval = 1000; // Store snapshot of pbwt
  int32_t nsamples=0, M=0;
  int32_t chunksize = 1000000;
  int32_t min_ac = 10;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("chr",&chrom, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("chunk-size",&chunksize, "Output by chunk")
    LONG_INT_PARAM("min-ac",&min_ac, "Minimum allele count to include")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file prefix")
    LONG_INT_PARAM("store-per-markers",&store_interval, "Snapshot per x markers")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // bcf reader
  std::vector<GenomeInterval> intervals;
  std::vector<std::string> v;
  if ( !reg.empty() ) {
    parse_intervals(intervals, "", reg);
    split(v, ":-", reg);
    chrom = v[0];
  }
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();
  nsamples = bcf_hdr_nsamples(odr.hdr);
  M = nsamples * 2;

  // Initialize pbwt Cursor
  pbwtCursor pc(M, 1);

  // Output
  std::string dout = "", aout = "";

  int32_t *p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t *info_ac = NULL;
  int32_t  n_ac = 0;
  bool *y;
  y = new bool[M];

  notice("Started Reading VCF, identifying %d samples. Will store snapshot every %d markers.", nsamples, store_interval);

  int32_t ck = 0, k_in_ck = 0;
  std::string sub = chrom + ":" + std::to_string(ck*chunksize+1) + "-" + std::to_string((ck+1)*chunksize);
  std::string d_outf = out + "_" + sub + "_prefix.dmat";
  std::string a_outf = out + "_" + sub + "_prefix.amat";
  std::ofstream d_wf, a_wf;
  d_wf.open(d_outf);
  a_wf.open(a_outf);
  // read marker
  for (int32_t k=0; odr.read(iv); ++k) {
    if (!bcf_is_snp(iv)) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (info_ac[0] < min_ac || info_ac[0] > 2*nsamples - min_ac) {continue;}
    if (bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }

    if (iv->pos+1 > (ck+1) * chunksize) { // Set up a new output file
      while (iv->pos+1 > (ck+2) * chunksize) { // Centromere
        ck++;
      }
      d_wf.close();
      a_wf.close();
      ck++;
      k_in_ck = 0;
      sub = chrom + ":" + std::to_string(ck*chunksize+1) + "-" + std::to_string((ck+1)*chunksize);
      d_outf = out + "_" + sub + "_prefix.dmat";
      a_outf = out + "_" + sub + "_prefix.amat";
      d_wf.open(d_outf);
      a_wf.open(a_outf);
    }

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
    pc.ForwardsAD_prefix(y, iv->pos+1);
    if (k_in_ck == 0 || k_in_ck % store_interval == 0) {

      d_wf << chrom << '\t' << iv->pos+1;
      for (int32_t j = 0; j < M; ++j)
        d_wf << '\t' << pc.d[j];
      d_wf << '\n';

      a_wf << chrom << '\t' << iv->pos+1;
      for (int32_t j = 0; j < M; ++j)
        a_wf << '\t' << pc.a[j];
      a_wf << '\n';
    }

    k_in_ck ++;

    if ( k % verbose == 0 )
    notice("Processed %d markers at chunk %d: %s:%d; snapshot at %d positions.", k, ck, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, k/store_interval);

  }
  d_wf.close();
  a_wf.close();

  delete [] y;

  return 0;
}



