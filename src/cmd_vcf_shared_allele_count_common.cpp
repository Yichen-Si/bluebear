#include "utils.h"
#include "cramore.h"
#include "hts_utils.h"

#include "compact_matrix.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"

#include <Eigen/Sparse>

int32_t cmdVcfXXtCommon(int32_t argc, char** argv) {

  std::string inVcf, reg, out, samples;
  int32_t verbose = 10000;
  bool haploid = false;
  bool count_alt = false;
  bool snp_only = false;
  bool annotate_carrier = false;
  double minAF = -1, maxAF = -1;
  int32_t minAC = -1, maxAC = -1, anno_max_ac = 10;

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
    LONG_PARAM("snp-only",&snp_only,"Ignore indels")
    LONG_DOUBLE_PARAM("min-af",&minAF,"Minimum minor/alt allele frequency")
    LONG_DOUBLE_PARAM("max-af",&maxAF,"Maximum minor/alt allele frequency")
    LONG_INT_PARAM("min-ac",&minAC,"Minimum minor/alt allele count")
    LONG_INT_PARAM("max-ac",&maxAC,"Maximum minor/alt allele count")
    LONG_PARAM("haploid",&haploid,"Input is haploid though vcf is read as diploid (?)")

	LONG_PARAM_GROUP("Output Options", NULL)
	LONG_STRING_PARAM("out", &out, "Output file prefix")
	LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")
    LONG_PARAM("count-alt",&count_alt,"Count the shared alt allele instead of minor allele (default)")
    LONG_PARAM("anno-carrier",&annotate_carrier,"Annotate carrier ID in an output vcf site file. Only for rare alleles")
    LONG_INT_PARAM("anno-max-ac",&anno_max_ac,"")


  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() ) {
	error("[E:%s:%d %s] --in-vcf, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

    std::vector<GenomeInterval> intervals;
    if ( !reg.empty() ) {
        parse_intervals(intervals, "", reg);
    }

    BCFOrderedReader odr(inVcf, intervals);
    bcf1_t* iv = bcf_init();

    if (!samples.empty()) {
        bcf_hdr_set_samples(odr.hdr, samples.c_str(), 1);
    }
    int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
    notice("Detected %d samples.", nsamples);

    int32_t nalleles = 2*nsamples;
    if (haploid) {nalleles = nsamples;}
    if (minAF >= 0) {
        minAC = (int32_t) minAF * nalleles;
    }
    if (maxAF >= 0 && maxAF <= 0.5) {
        maxAC = (int32_t) maxAF * nalleles;
    }
    if (minAC < 0) { minAC = 1; }
    if (maxAC < 0) { maxAC = (int32_t) (nalleles / 2); }
    if (anno_max_ac < minAC) {annotate_carrier = false; }

    std::string outf = out + ".site.vcf.gz";
    BCFOrderedWriter odw(outf.c_str(),0);
    bcf_hdr_t* hnull = bcf_hdr_subset(odr.hdr, 0, 0, 0);
    bcf_hdr_remove(hnull, BCF_HL_FMT, NULL);
    odw.set_hdr(hnull);
    char buffer[65536];
	if ( annotate_carrier && bcf_hdr_id2int(odw.hdr, BCF_DT_ID, "CarrierID") < 0 ) {
		sprintf(buffer,"##INFO=<ID=CarrierID,Number=1,Type=String,Description=\"Minor allele carriers\">\n");
		bcf_hdr_append(odw.hdr, buffer);
	}
    odw.write_hdr();

    // handle filter string
    vfilt.init(odr.hdr);
    if ( snp_only ) vfilt.snpOnly = true;

    std::vector<int32_t> freq;
    bitmatrix bmat(nsamples);

    int32_t n_gt = 0, n_ac = 0, n_an = 0;
    int32_t* p_gt = NULL;
    int32_t *info_ac = NULL;
    int32_t *info_an = NULL;
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

        if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
        }
        if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
        if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
        int32_t ac = info_ac[0];
        if (info_an[0] > nalleles && haploid) {ac = ac / 2;}
        if (ac > nalleles / 2) {ac = nalleles - ac;}
        if (ac < minAC || ac > maxAC) {continue;}

        ++nVariant;

        bool flip = 0;
        ac = 0;
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
        if ((!count_alt) && ac > nalleles / 2) {
            flip = 1;
            ac = nalleles - ac;
            for(int32_t i = 0; i < nalleles; ++i) {
                gt[i] = 1 - gt[i];
            }
        }
        // freq.push_back((double) ac/nalleles);
        freq.push_back(ac);
        bmat.add_row_bytes(gt);


        bcf1_t* nv = bcf_dup(iv);
        bcf_unpack(nv, BCF_UN_ALL);
        bcf_subset(odw.hdr, nv, 0, 0);
        odw.write(nv);
        bcf_destroy(nv);

        // char** alleles = bcf_get_allele(iv);
        // hprintf(vwf, "%s\t%d\t%s\t%s\t%d\t%d\n", bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, alleles[0], alleles[1], (int32_t) flip, ac);
    }
    odr.close();
    odw.close();
    // hts_close(vwf);

    notice("Matrix: %d x %d.\n", bmat.nrow, bmat.ncol);

    outf = out + ".count.mtx.gz";
    htsFile* wf = hts_open(outf.c_str(), "wg");
    for (int32_t i = 0; i < bmat.nrow - 1; ++i) {
        bool rare = 0;
        if (freq[i] < output_oneside_minAC) {
            rare = 1;
        }
        for (int32_t j = i; j < bmat.nrow; ++j) {
            if (rare && freq[j] < output_oneside_minAC) {
                continue;
            }
            int32_t mmt = bmat.inner_prod_and_bytes(i, j);
            if (mmt < 1) {continue;}
            if (count_only) {
                hprintf(wf, "%d\t%d\t%d\n", i, j, mmt);
            } else {
                double xy = (double) mmt / nalleles;
                double f1 = (double) freq[i] / nalleles;
                double f2 = (double) freq[j] / nalleles;
                double d  = xy - f1 * f2;
                double dmax = 1;
                if (d < 0) {
                    dmax = std::max( -f1*f2, -(1-f1)*(1-f2) ) ;
                }
                if (d > 0) {
                    dmax = std::min( f1*(1-f2), (1-f1)*f2 ) ;
                }
                double rsq = d * d / f2 / (1-f1) / f1 / (1-f2);
                d = d / dmax;
                if (d < minDprime) {
                    continue;
                }
                hprintf(wf, "%d\t%d\t%d\t%.4f\t%.4f\n", i, j, mmt, d, rsq);
            }
        }
        if (i % 100 == 0)
            notice("Outputing pairs... one side %d / %d.\n", i, bmat.nrow);
    }

    hts_close(wf);

    return 0;
}
