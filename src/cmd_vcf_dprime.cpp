#include "utils.h"
#include "cramore.h"
#include "hts_utils.h"

#include "compact_matrix.h"
#include "bcf_filter_arg.h"
// #include "bcf_filtered_reader.h"
#include "bcf_ordered_reader.h"

int32_t cmdVcfDprime(int32_t argc, char** argv) {

  std::string inVcf, reg, out, samples;
  int32_t verbose = 10000;
  double minDprime = -1;
  bool haploid = false;
  bool count_only = false;
  bool recode_minor = false;
  bool includeIndels = false;

  bcf_vfilter_arg vfilt;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
	LONG_PARAM_GROUP("Input Sites", NULL)
	LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
	LONG_STRING_PARAM("samples",&samples,"List of sample ID")

	LONG_PARAM_GROUP("Variant Filtering Options", NULL)
	LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
	LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
	LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")
    LONG_PARAM("include-indels",&includeIndels,"Include indels")

	LONG_PARAM_GROUP("Output Options", NULL)
	LONG_STRING_PARAM("out", &out, "Output file prefix")
	LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")
    LONG_DOUBLE_PARAM("minDprime",&minDprime,"Minimum Dprime of output pairs")
    LONG_PARAM("haploid",&haploid,"Input is haploid though vcf is read as diploid (?)")
    LONG_PARAM("count-only",&count_only,"If only compute genotype inner product without calculating statistics")
    LONG_PARAM("recode-minor",&recode_minor,"If the number of overlapped carriers is always calculated for the minor alleles - cannot be used when only partial data is available")


  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() ) {
	error("[E:%s:%d %s] --in-vcf, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

    std::string outf = out + ".var.tsv";
    htsFile* vwf = hts_open(outf.c_str(), "w");
    outf = out + ".dprime.mtx.gz";
    htsFile* dwf = hts_open(outf.c_str(), "wg");

    std::vector<GenomeInterval> intervals;
    if ( !reg.empty() ) {
        parse_intervals(intervals, "", reg);
    }

    // BCFFilteredReader odr(inVcf);
    BCFOrderedReader odr(inVcf, intervals);
    bcf1_t* iv = bcf_init();

    if (!samples.empty()) {
        bcf_hdr_set_samples(odr.hdr, samples.c_str(), 1);
    }

    // handle filter string
    vfilt.init(odr.hdr);
    if ( !includeIndels ) vfilt.snpOnly = true;
    int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
    notice("Processing %d samples.", nsamples);
    int32_t nalleles = 2*nsamples;
    if (haploid) {nalleles = nsamples;}

    std::vector<double> freq;
    bitmatrix bmat(nsamples);

    int32_t* p_gt = NULL;
    int32_t n_gt = 0;
    int32_t nVariant = 0;

    for(int32_t k=0; odr.read(iv); ++k) {  // read marker

        if ( k % verbose == 0 )
            notice("Processing %d markers at %s:%d.", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);

        bcf_unpack(iv, BCF_UN_FLT);
        if ( vfilt.snpOnly && (!bcf_is_snp(iv)) ) {continue;}


        // check --apply-filters
        bool has_filter = vfilt.req_flt_ids.empty() ? true : false;
        if ( ! has_filter ) {
            for(int32_t i=0; i < iv->d.n_flt; ++i) {
                for(int32_t j=0; j < (int32_t)vfilt.req_flt_ids.size(); ++j) {
                    if ( vfilt.req_flt_ids[j] == iv->d.flt[i] )
                        has_filter = true;
                }
            }
        }
        if ( ! has_filter ) { continue; }
        // check filter logic
        if ( vfilt.filt != NULL ) {
            int32_t ret = filter_test(vfilt.filt, iv, NULL);
            if ( vfilt.filter_logic == FLT_INCLUDE ) { if ( !ret)  has_filter = false; }
            else if ( ret ) { has_filter = false; }
        }
        if ( ! has_filter ) { continue; }

        // if (!ord.passed_vfilter()) {
        //     continue;
        // }

        if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
        }
        ++nVariant;

        bool flip = 0;
        int32_t ac = 0;
        uint8_t gt[nalleles];
        if (haploid) {
            for(int32_t i=0; i < nsamples; ++i) {
                uint8_t geno = 0;
                if ( bcf_gt_is_missing(p_gt[i*2]) ) {
                    geno = 0;
                }
                else {
                    geno = (bcf_gt_allele(p_gt[i*2]) > 0) ? 1 : 0;
                }
                gt[i] = geno;
                ac += geno;
            }
        } else {
            for(int32_t i=0; i < n_gt; ++i) {
                uint8_t geno = 0;
                if ( bcf_gt_is_missing(p_gt[i]) ) {
                    geno = 0;
                }
                else {
                    geno = (bcf_gt_allele(p_gt[i]) > 0) ? 1 : 0;
                }
                gt[i] = geno;
                ac += geno;
            }
        }
        if (recode_minor && ac > nalleles / 2) {
            flip = 1;
            ac = nalleles - ac;
            for(int32_t i = 0; i < nalleles; ++i) {
                gt[i] = 1 - gt[i];
            }
        }
        freq.push_back((double) ac/nalleles);
        bmat.add_row_bytes(gt);
        char** alleles = bcf_get_allele(iv);
        hprintf(vwf, "%s\t%d\t%s\t%s\t%d\t%d\n", bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, alleles[0], alleles[1], (int32_t) flip, ac);
    }
    odr.close();

    notice("Matrix: %d x %d.\n", bmat.nrow, bmat.ncol);

    for (int32_t i = 0; i < bmat.nrow - 1; ++i) {
        for (int32_t j = i; j < bmat.nrow; ++j) {
            int32_t mmt = bmat.inner_prod_and_bytes(i, j);
            if (mmt < 1) {continue;}

            if (count_only) {
                hprintf(dwf, "%d\t%d\t%d\n", i, j, mmt);
            } else {
                double xy = (double) mmt / nalleles;
                double d  = xy - freq[i] * freq[j];
                double dmax = 1;
                if (d < 0) {
                    dmax = std::max( -freq[i]*freq[j], -(1-freq[i])*(1-freq[j]) ) ;
                }
                if (d > 0) {
                    dmax = std::min( freq[i]*(1-freq[j]), (1-freq[i])*freq[j] ) ;
                }
                double rsq = d * d / freq[j] / (1-freq[i]) / freq[i] / (1-freq[j]);
                d = d / dmax;
                if (d < minDprime) {
                    continue;
                }
                hprintf(dwf, "%d\t%d\t%d\t%.4f\t%.4f\n", i, j, mmt, d, rsq);
            }
        }
    }

    hts_close(dwf);
    hts_close(vwf);

    return 0;
}
