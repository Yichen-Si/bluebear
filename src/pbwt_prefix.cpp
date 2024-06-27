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
// TODO: handle missing genotypes better
int32_t pbwtBuildPrefix(int32_t argc, char** argv) {

  std::string inVcf, chrom, reg, out;
  int32_t verbose = 10000;
  int32_t max_store_interval = 100000; // bp
  int32_t store_interval = 10000;      // markers
  int32_t min_ac = 10;
  double max_missing = 0.05;
  bool haploid = false;
  bool snp_only = false;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("chr",&chrom, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("min-ac",&min_ac, "Minimum allele count to include")
		LONG_PARAM("haploid",&haploid, "Assume input is haploid coded in diploid vcf form and use only the first allele e.g. for chrY")
		LONG_PARAM("snp-only",&snp_only, "Only consider SNPs")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file prefix")
    LONG_INT_PARAM("store-per-markers",&store_interval, "Snapshot per X markers (after filtering)")
    LONG_INT_PARAM("max-store-interval",&max_store_interval, "Snapshot at least per X bp")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();
	if ( inVcf.empty() || out.empty() ) {
		error("[E:%s:%d %s] --in-vcf, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
	}

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
  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
  int32_t M = nsamples;
  if (!haploid) {
    M *= 2;
  }
	if (nsamples == 0) {
		error("Did not find any samples in the VCF file.");
	}
  notice("Identifying %d samples. Will store snapshot every %d markers or after %d bp", nsamples, store_interval, max_store_interval);

  // Initialize pbwt Cursor
  pbwtCursor pc(M, 1);

  // Output
  std::string d_outf = out + ".prefix.dmat";
  std::string a_outf = out + ".prefix.amat";
  std::ofstream d_wf, a_wf;
  d_wf.open(d_outf);
  a_wf.open(a_outf);

  int32_t *p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t *info_ac = NULL, *info_an = NULL;
  int32_t  n_ac = 0, n_an = 0;
  bool* y = new bool[M];
  int32_t last_stored_pos = 0;

  // read marker
  int32_t k = 0, n_rec = 0, n_out = 0;
  for (k=0; odr.read(iv); ++k) {
    if ( k % verbose == 0 ) {
      notice("Processed %d markers %s:%d; used %d; snapshot at %d positions.", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, n_rec, n_out);
    }
    if (snp_only && !bcf_is_snp(iv)) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
    if (info_ac[0] < min_ac || info_ac[0] > info_an[0] - min_ac) {continue;}
    if (bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }
    if (chrom.empty()) {
      chrom = std::string(bcf_hdr_id2name(odr.hdr, iv->rid));
    }

    int32_t major = (info_ac[0] > info_an[0]/2) ? 1 : 0;
    int32_t n_miss = 0;
    for (int32_t i = 0; i < M; ++i) {
      int32_t gt = p_gt[i];
      if (haploid) {
        gt = p_gt[i*2];
      }
      if (bcf_gt_is_missing(gt)) {
        n_miss++;
        y[i] = major;
      } else {
        y[i] = ((bcf_gt_allele(gt) > 0) ? 1 : 0);
      }
    }
    if (n_miss > max_missing * M) {continue;}

    if (n_rec == 0) { // proper init given the first *observed* position
      pc.Reset(M, iv->pos+1);
    }
    pc.ForwardsAD_prefix(y, iv->pos+1);
    if (n_rec % store_interval == 0 || iv->pos - last_stored_pos > max_store_interval) {
      d_wf << chrom << '\t' << iv->pos+1;
      for (int32_t j = 0; j < M; ++j)
        d_wf << '\t' << pc.d[j];
      d_wf << '\n';
      a_wf << chrom << '\t' << iv->pos+1;
      for (int32_t j = 0; j < M; ++j)
        a_wf << '\t' << pc.a[j];
      a_wf << '\n';
      last_stored_pos = iv->pos;
      n_out++;
    }
    n_rec++;
  }
  d_wf.close();
  a_wf.close();

  delete [] y;
  delete [] p_gt;
  bcf_destroy(iv);
  odr.close();

  notice("Finished. Processed %d markers, used %d, snapshot at %d positions.", k, n_rec, n_out);
  return 0;
}
