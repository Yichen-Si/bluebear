#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "compact_matrix.h"

#include <functional>
#include <iomanip>


int32_t test(int32_t argc, char** argv) {
  std::string inVcf, out, reg;
  std::string outf;
  std::string chrom="chr20";
  int32_t verbose = 1000;
  int32_t nsamples=0, M=0;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
    LONG_STRING_PARAM("chr",&chrom,"Chromosome to process")

    LONG_PARAM_GROUP("Additional Options", NULL)

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file prefix")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

    notice("Processing %s.", reg.c_str());
    std::vector<GenomeInterval> intervals;
    parse_intervals(intervals, "", reg);
    BCFOrderedReader odr(inVcf, intervals);
    nsamples = bcf_hdr_nsamples(odr.hdr);
    M = nsamples * 2;
    bcf1_t* iv = bcf_init();

    outf = out + ".flipped.bcf";
    htsFile *fp = NULL;

    char mode[] = "wb";
    fp = hts_open(outf.c_str(), &mode[0]);
    bcf_hdr_write(fp, odr.hdr);

    int32_t* p_gt = NULL;
    int32_t n_gt = 0;

    for (int32_t k=0; odr.read(iv); ++k) { // Read haplotypes
      if (bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
      }
      int32_t y[M];
      for (int32_t i = 0; i <  nsamples; ++i) {
        int32_t g1 = p_gt[2*i];
        int32_t g2 = p_gt[2*i+1];
        if (bcf_gt_is_missing(g1)) {
          y[2*i] = 0;
        } else {
          y[2*i] = ((bcf_gt_allele(g1) > 0) ? 0:1) ;
        }
        if (bcf_gt_is_missing(g2)) {
          y[2*i+1] = 0;
        } else {
          y[2*i+1] = 1;
        }
        y[2*i]   = bcf_gt_phased(y[2*i]);
        y[2*i+1] = bcf_gt_phased(y[2*i+1]);
        // std::cout << bcf_gt_allele(g1) << "|" << bcf_gt_allele(g2) << "->" << y[2*i] << "|" << y[2*i+1] << ' ';
      } // Finish reading one SNP
      // std::cout << '\n';
      if (bcf_update_genotypes(odr.hdr, iv, y, M)) {
        error("Cannot update GT\n");
      }
      bcf_write(fp, odr.hdr, iv);
    } // Finish reading all haplotypes in this chunk

  hts_close(fp);

  return 0;
}














