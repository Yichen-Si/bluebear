#include "utils.h"
#include "cramore.h"
#include "hts_utils.h"

#include "compact_matrix.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"

#include <Eigen/Sparse>

int32_t cmdVcfXXtRare(int32_t argc, char** argv) {

  std::string inVcf, reg, out, samples;
  int32_t verbose = 10000;
  bool haploid = false;
  bool count_alt = false;
  bool snp_only = false;
  bool annotate_carrier = false;
  double minAF = -1, maxAF = -1;
  int32_t minAC = -1, maxAC = -1, anno_max_ac = 10;
  std::string tag = "";

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
    LONG_STRING_PARAM("key-of-weight",&tag, "INFO/TAG as multiples/weights to add up")
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
    Eigen::SparseMatrix<int32_t > XXt(nsamples, nsamples);
    XXt.reserve(Eigen::VectorXi::Constant(nsamples,maxAC*5));

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
		sprintf(buffer,"##INFO=<ID=CarrierID,Number=.,Type=Integer,Description=\"Indecies of minor allele carriers\">\n");
		bcf_hdr_append(odw.hdr, buffer);
	}
    odw.write_hdr();

    bool use_weight = false;
    if ( bcf_hdr_id2int(odw.hdr, BCF_DT_ID, tag.c_str()) >= 0 ) {
        use_weight = true;
        notice("Detected TAG %s to use as weights of shared alleles", tag.c_str());
    }

    // handle filter string
    vfilt.init(odr.hdr);
    if ( snp_only ) vfilt.snpOnly = true;

    int32_t n_gt = 0, n_ac = 0, n_an = 0, n_mult = 0;
    int32_t *p_gt = NULL;
    int32_t *info_ac = NULL, *info_an = NULL, *info_mult = NULL;
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
        int32_t weight = 1;
        if (use_weight && bcf_get_info_int32(odr.hdr, iv, tag.c_str(), &info_mult, &n_mult) >= 0) {
            weight = info_mult[0];
        }

        int32_t ac = info_ac[0];
        if (info_an[0] > nalleles && haploid) {ac = ac / 2;}
        int32_t minor = 1;
        if (ac > nalleles / 2) {
            ac = nalleles - ac;
            minor = 0;
        }
        if (ac < minAC || ac > maxAC) {continue;}

        ++nVariant;
        std::vector<int32_t> carry;
        for(int32_t i=0; i < nsamples; ++i) {
            uint8_t geno = 0;
            if (haploid) {
                if ( bcf_gt_is_missing(p_gt[i*2]) )
                    continue;
                geno = (bcf_gt_allele(p_gt[i*2]) > 0) ? 1 : 0;
                geno = (geno == minor) ? 1 : 0;
            } else {
                if ( bcf_gt_is_missing(p_gt[i*2]) || bcf_gt_is_missing(p_gt[i*2+1])  ) {
                    continue;
                }
                if (minor == 1) {
                    geno = ((bcf_gt_allele(p_gt[i*2]) > 0 || bcf_gt_allele(p_gt[i*2+1]) > 0) ? 1 : 0);
                } else {
                    geno = ((bcf_gt_allele(p_gt[i*2]) == 0 || bcf_gt_allele(p_gt[i*2+1]) == 0) ? 1 : 0);
                }
            }
            if (geno) {
                carry.push_back(i);
            }
        }
        uint32_t ncarry = carry.size();
        for (uint32_t i = 0; i < ncarry; i ++) {
            for (uint32_t j = i; j < ncarry; j++) {
                XXt.coeffRef(carry[i],carry[j]) += weight;
            }
        }
        bcf1_t* nv = bcf_dup(iv);
        bcf_unpack(nv, BCF_UN_ALL);
        bcf_subset(odw.hdr, nv, 0, 0);
        if (ac <= anno_max_ac) {
            bcf_update_info_int32(odw.hdr, nv, "CarrierID", &(carry[0]), ncarry);
        }
        odw.write(nv);
        bcf_destroy(nv);
    }
    odr.close();
    odw.close();
    XXt.makeCompressed();
    notice("Recorded %d non-zero elements including diagonals.\n", XXt.nonZeros());

    outf = out + ".XXt.gz";
    htsFile* wf = hts_open(outf.c_str(), "wg");

    for (int32_t i = 0; i < nsamples; ++i) {
        for (Eigen::SparseMatrix<int32_t>::InnerIterator it(XXt,i); it; ++it) {
            hprintf(wf, "%d\t%d\t%d\n", it.row(), it.col(), it.value());
        }
    }
    hts_close(wf);

    return 0;
}
