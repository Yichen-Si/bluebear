#include "utils.h"
#include "cramore.h"
#include "hts_utils.h"

#include "compact_matrix.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"

#include "binary_variant_compress.h"

double dprime(int32_t n, int32_t n1x, int32_t nx1, int32_t n11) {
    double p1  = (double)n1x / n;
    double p2  = (double)nx1 / n;
    double p12 = (double)n11 / n;
    double d = p12 - p1 * p2;
    double dmax = 0;
    if (d < 0) {
        dmax = std::min(p1 * p2, (1 - p1) * (1 - p2));
    } else {
        dmax = std::min(p1 * (1 - p2), (1 - p1) * p2);
    }
    return d / dmax;
}

int32_t cmdDprimeFromBinaryStorage(int32_t argc, char** argv) {

    std::string inBin, out, outf;
    int32_t verbose = 10000;
    int32_t debug = 0;
    int32_t acmin = 2;

    paramList pl;
    BEGIN_LONG_PARAMS(longParameters)
        LONG_PARAM_GROUP("Input Sites", NULL)
        LONG_STRING_PARAM("in-bin",&inBin, "Input binary file")
        LONG_PARAM_GROUP("Output Options", NULL)
        LONG_STRING_PARAM("out", &out, "Output file prefix")
        LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")
        LONG_INT_PARAM("debug",&debug,"Debug")
        LONG_INT_PARAM("min-ac",&acmin,"Minimum allele count")
    END_LONG_PARAMS();
    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    std::vector<BinaryVariant> mixStorage;
    parseBinaryFile(inBin, mixStorage);
    int32_t nValues = mixStorage.size();
    int32_t nalleles = mixStorage[0].nBits;
    notice("Loaded %d variants", nValues);

    outf = out + ".dprime.pairs.tsv.gz";
    htsFile* wf = hts_open(outf.c_str(), "wz");
    if (wf == NULL) {
        error("Could not open file: %s", outf.c_str());
    }

    int32_t npairs = 0, nrec = 0;
    for (int32_t i = 1; i < nValues; ++i) {
        if (mixStorage[i].mOnes < acmin) {
            npairs += i;
            continue;
        }
        int32_t n1x = mixStorage[i].mOnes;
        for (int32_t j = 0; j < i; ++j) {
            npairs++;
            if (mixStorage[j].mOnes < acmin) {
                continue;
            }
            int32_t nx1 = mixStorage[j].mOnes;
            int32_t n11 = binary_intersect(mixStorage[i], mixStorage[j]);
            double dp = dprime(nalleles, n1x, nx1, n11);
            hprintf(wf, "%d\t%d\t%d\t%d\t%d\t%.6f\n", i, j, n1x, nx1, n11, dp);
            nrec++;
            if (nrec % verbose == 0) {
                notice("Processed %d pairs, recorded %d", npairs, nrec);
            }
        }
        if (debug && nrec > debug) {
            break;
        }
    }

    notice("Processed %d pairs, recorded %d", npairs, nrec);
    hts_close(wf);

    return 0;

}














int32_t cmdVcfDprime(int32_t argc, char** argv) {

  std::string inVcf, reg, out, samples;
  int32_t verbose = 10000;
  double minDprime = -1;
  bool debug = false;
  bool haploid = false;
  bool count_only = false;
  bool recode_minor = false;
  bool all_pair_fail_fourgamete = false;
  bool includeIndels = false;
  double minAF = -1, maxAF = -1;
  int32_t minAC = -1, maxAC = -1;
  int32_t output_oneside_minAC = -1;

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
    LONG_DOUBLE_PARAM("min-af",&minAF,"Minimum minor allele frequency")
    LONG_DOUBLE_PARAM("max-af",&maxAF,"Maximum minor allele frequency")
    LONG_INT_PARAM("min-ac",&minAC,"Minimum minor allele count")
    LONG_INT_PARAM("max-ac",&maxAC,"Maximum minor allele count")
    LONG_INT_PARAM("out-oneside-min-ac",&output_oneside_minAC,"Only output pairs with at least one variant with minor allele count larger than")
    LONG_PARAM("haploid",&haploid,"Input is haploid though vcf is read as diploid (?)")

	LONG_PARAM_GROUP("Output Options", NULL)
	LONG_STRING_PARAM("out", &out, "Output file prefix")
	LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")
    LONG_DOUBLE_PARAM("minDprime",&minDprime,"Minimum Dprime of output pairs")
    LONG_PARAM("count-only",&count_only,"If only compute genotype inner product without calculating statistics")
    LONG_PARAM("recode-minor",&recode_minor,"If the number of overlapped carriers is always calculated for the minor alleles - cannot be used when only partial data is available")
    LONG_PARAM("four-gamete",&all_pair_fail_fourgamete,"Output all and only pairs that fail four-gamete test")

    LONG_INT_PARAM("debug",&debug,"")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() ) {
	error("[E:%s:%d %s] --in-vcf, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

    std::string outf = out + ".dprime.mtx.gz";
    htsFile* wf = hts_open(outf.c_str(), "wg");

    std::vector<GenomeInterval> intervals;
    if ( !reg.empty() ) {
        parse_intervals(intervals, "", reg);
    }

    BCFOrderedReader odr(inVcf, intervals);
    bcf1_t* iv = bcf_init();

    int32_t nsamples_full = bcf_hdr_nsamples(odr.hdr);
    if (!samples.empty()) {
        bcf_hdr_set_samples(odr.hdr, samples.c_str(), 1);
    }
    int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
    notice("Detected %d samples.", nsamples);

    outf = out + ".site.vcf.gz";
    BCFOrderedWriter odw(outf.c_str(),0);
    bcf_hdr_t* hnull = bcf_hdr_subset(odr.hdr, 0, 0, 0);
    bcf_hdr_remove(hnull, BCF_HL_FMT, NULL);
    odw.set_hdr(hnull);
    odw.write_hdr();

    // handle filter string
    vfilt.init(odr.hdr);
    if ( !includeIndels ) vfilt.snpOnly = true;
    int32_t nalleles = 2*nsamples;
    int32_t nalleles_full = 2*nsamples_full;
    if (haploid) {
        nalleles = nsamples;
        nalleles_full = nsamples_full;
    }
    if (minAF >= 0) {
        minAC = (int32_t) (minAF * nalleles_full);
    }
    if (maxAF >= 0 && maxAF <= 0.5) {
        maxAC = (int32_t) (maxAF * nalleles_full);
    }
    if (minAC < 0) { minAC = 1; }
    if (maxAC < 0) { maxAC = (int32_t) (nalleles_full / 2); }
    if (output_oneside_minAC < 0) {output_oneside_minAC = minAC; }
    if (debug) {
        std::cout << minAC << '\t' << maxAC << '\n';
    }

    std::vector<int32_t> freq;
    bitmatrix bmat(nalleles);

    int32_t n_gt = 0, n_ac = 0, n_an = 0;
    int32_t* p_gt = NULL;
    int32_t *info_ac = NULL;
    int32_t *info_an = NULL;
    int32_t nVariant = 0;

    for(int32_t k=0; odr.read(iv); ++k) {  // read marker

        if ( k % verbose == 0 )
            notice("Processed %d markers at %s:%d.", nVariant, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);

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
        // Filter based on original AC/AN in case of sub-sampling
        int32_t minor_ac = info_ac[0];
        bool flip = false;
        if (info_an[0] > nalleles_full && haploid) {minor_ac = minor_ac / 2;}
        if (minor_ac > info_an[0]/2) {
            flip = true;
            minor_ac = info_an[0] - minor_ac;
        }
        if (minor_ac < minAC || minor_ac > maxAC) {continue;}
        ++nVariant;
        if (debug) {
            std::cout << info_an[0] << '\t' << minor_ac << '\t' << minor_ac*1./info_an[0] << '\n';
        }

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
            }
        }
        if (recode_minor && flip) {
            for(int32_t i = 0; i < nalleles; ++i) {
                gt[i] = 1 - gt[i];
            }
        }
        freq.push_back(minor_ac);
        bmat.add_row_bytes(gt);

        bcf1_t* nv = bcf_dup(iv);
        bcf_unpack(nv, BCF_UN_ALL);
        bcf_subset(odw.hdr, nv, 0, 0);
        odw.write(nv);
        bcf_destroy(nv);

        // if (debug && nVariant > 100) {exit(0);}
    }
    odr.close();
    odw.close();

    notice("Matrix: %d x %d.\n", bmat.nrow, bmat.ncol);

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
            if (mmt < 1) {continue;} // Missing 11
            if (count_only) {
                hprintf(wf, "%d\t%d\t%d\n", i, j, mmt);
            } else {
                if (all_pair_fail_fourgamete) {
                    if (freq[i] == mmt || freq[j] == mmt) {
                        continue; // Missing 10 or 01
                    }
                    if (freq[i]+freq[j]-mmt == nalleles) {
                        continue; // Missing 00
                    }
                }
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
