#include "cramore.h"
#include "hts_utils.h"
#include "utils.h"
#include <algorithm>
#include <iostream>
#include <cmath>
#include "fstream"
#include "sstream"

// #include <Eigen/Dense>
// using Eigen::MatrixXd;
// using Eigen::VectorXd;
// using Eigen::ArrayXXd;

// #include "fa_reader.h"
// #include "seq_basics.h"

#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"

int32_t test(int32_t argc, char** argv) {

  std::string inVcf, reg, out, samples;
  std::vector<std::string> test_samples, ctrl_samples;
  std::vector<uint32_t> test_sample_idx, ctrl_sample_idx;
  bool preview = false;
  int32_t verbose = 10000;
  int32_t min_ctrl_support = 1;
  int32_t max_minor_adp_ctrl = 2;
  double  max_minor_vaf_ctrl = .1;
  double  min_minor_vaf_test = .2;
  int32_t min_minor_adp_test = 5;
  int32_t max_sum_minor_abq_ctrl = 60;
  int32_t minDP = 8;
  double  tol = 1e-8;
  int32_t max_q = 80;
  bool single_strand = false; // if true, only use reads from one strand

  bcf_vfilter_arg vfilt;

  paramList pl;
  BEGIN_LONG_PARAMS(longParameters)
	LONG_PARAM_GROUP("Input Information", NULL)
	LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
	LONG_STRING_PARAM("samples",&samples,"List of sample ID")
    LONG_MULTI_STRING_PARAM("test-samples",&test_samples,"List of test sample ID")
    LONG_MULTI_STRING_PARAM("ctrl-samples",&ctrl_samples,"List of reference sample ID")
    LONG_MULTI_INT_PARAM("test-sample-index",&test_sample_idx,"List of test sample indices")
    LONG_MULTI_INT_PARAM("ctrl-sample-index",&ctrl_sample_idx,"List of reference sample indices")
    LONG_PARAM("single-strand",&single_strand,"Use only reads from one strand")

	LONG_PARAM_GROUP("Variant Filtering Options", NULL)
	LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")
    LONG_INT_PARAM("minDP",&minDP,"Min DP for each sample")
    LONG_INT_PARAM("min-ctrl-support",&min_ctrl_support,"Min number of reference samples with called homozygouse genotype to declare a potential test-specific allele")
    LONG_INT_PARAM("max-minor-adp-ctrl",&max_minor_adp_ctrl,"")
    LONG_DOUBLE_PARAM("max-minor-vaf-ctrl",&max_minor_vaf_ctrl,"")
    LONG_INT_PARAM("max-sum-minor-abq-ctrl",&max_sum_minor_abq_ctrl,"")

	LONG_PARAM_GROUP("Output Options", NULL)
	LONG_STRING_PARAM("out", &out, "Output file prefix")
	LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")
    LONG_PARAM("preview", &preview, "")
    LONG_DOUBLE_PARAM("min-minor-vaf-test",&min_minor_vaf_test,"For preview only")
    LONG_INT_PARAM("min-minor-adp-test",&min_minor_adp_test,"For preview only")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() ) {
	error("[E:%s:%d %s] --in-vcf, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }
  if ( test_samples.empty() && test_sample_idx.empty() ) {
    error("[E:%s:%d %s] --test-samples or --test-sample-index is required",__FILE__,__LINE__,__FUNCTION__);
  }

    std::string outf = out + ".tsv.gz";
    htsFile* wf = hts_open(outf.c_str(), "wg");

    std::vector<GenomeInterval> intervals;
    if ( !reg.empty() ) {
        parse_intervals(intervals, "", reg);
    }
    BCFOrderedReader odr(inVcf, intervals);
    bcf1_t* iv = bcf_init();
    if (!samples.empty()) {
        bcf_hdr_set_samples(odr.hdr, samples.c_str(), 1);
    }
    // Parse input sample groups
    int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
    if (test_sample_idx.empty() && !test_samples.empty()) {
        // Translate input sample ID to sample indices
        for (auto & v : test_samples) {
            int idx = bcf_hdr_id2int(odr.hdr,BCF_DT_SAMPLE,v.c_str());
            if (idx >= 0 && idx < nsamples) {
                test_sample_idx.push_back(idx);
            }
            else {
                warning("[W:%s:%d %s] Cannot find sample %s in the VCF header",__FILE__,__LINE__,__FUNCTION__,v.c_str());
            }
        }
    }
    if (ctrl_sample_idx.empty() && !ctrl_samples.empty()) {
        for (auto & v : ctrl_samples) {
            int idx = bcf_hdr_id2int(odr.hdr,BCF_DT_SAMPLE,v.c_str());
            if (idx >= 0 && idx < nsamples) {
                ctrl_sample_idx.push_back(idx);
            }
            else {
                warning("[W:%s:%d %s] Cannot find sample %s in the VCF header",__FILE__,__LINE__,__FUNCTION__,v.c_str());
            }
        }
    }
    std::sort( test_sample_idx.begin(), test_sample_idx.end() );
    test_sample_idx.erase( std::unique( test_sample_idx.begin(), test_sample_idx.end() ), test_sample_idx.end() );
    if (!ctrl_sample_idx.empty()) {
        std::sort( ctrl_sample_idx.begin(), ctrl_sample_idx.end() );
        ctrl_sample_idx.erase( std::unique( ctrl_sample_idx.begin(), ctrl_sample_idx.end() ), ctrl_sample_idx.end() );
    }
    if (!test_sample_idx.empty() && ctrl_sample_idx.empty()) {
        // Assume all other samples are reference samples
        uint32_t j = 0;
        for (auto & v : test_sample_idx) {
            while (j < v) {
                ctrl_sample_idx.push_back(j);
                ++j;
            }
            ++j;
        }
        while (j < nsamples) {
            ctrl_sample_idx.push_back(j);
            ++j;
        }
    }
    int32_t ntest = test_sample_idx.size();
    int32_t nctrl = ctrl_sample_idx.size();
    notice("Processing %d samples (%d test v.s. %d reference).", nsamples, ntest, nctrl);

    // handle filter string
    std::string filter_str;
    int32_t filter_logic = 0;
    if ( vfilt.include_expr.empty() ) {
        if ( vfilt.exclude_expr.empty() ) {
            // do nothing
        }
        else {
            filter_str = vfilt.exclude_expr;
            filter_logic |= FLT_EXCLUDE;
        }
    }
    else {
        if ( vfilt.exclude_expr.empty() ) {
            filter_str = vfilt.include_expr;
            filter_logic |= FLT_INCLUDE;
        }
        else {
            error("[E:%s:%d %s] Cannot use both --include-expr and --exclude-expr options",__FILE__,__LINE__,__FUNCTION__);
        }
    }
    filter_t* filt = NULL;
    if ( filter_logic != 0 )
        filter_init(odr.hdr, filter_str.c_str());
    // handle --apply-filtrs
    std::vector<int32_t> req_flt_ids;
    if ( !vfilt.required_filters.empty() ) {
        for(int32_t i=0; i < (int32_t)vfilt.required_filters.size(); ++i) {
            req_flt_ids.push_back(bcf_hdr_id2int(odr.hdr, BCF_DT_ID, vfilt.required_filters[i].c_str()));
        }
    }

    int32_t* v_dp = NULL;
    int32_t* v_rbq = NULL;
    int32_t* v_abq = NULL;
    int32_t* v_rdf = NULL;
    int32_t* v_adf = NULL;
    int32_t* v_rdr = NULL;
    int32_t* v_adr = NULL;
    int32_t n_fmt = 0;
    int32_t nVariant = 0;

    for(int32_t k=0; odr.read(iv); ++k) {  // read marker

        if ( k > 0 && k % verbose == 0 )
            notice("Processing %d markers at %s:%d. Kept %d variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nVariant);

        bcf_unpack(iv, BCF_UN_FLT);
        // check --apply-filters
        bool has_filter = req_flt_ids.empty() ? true : false;
        if ( ! has_filter ) {
            for(int32_t i=0; i < iv->d.n_flt; ++i) {
                for(int32_t j=0; j < (int32_t)req_flt_ids.size(); ++j) {
                    if ( req_flt_ids[j] == iv->d.flt[i] )
                        has_filter = true;
                }
            }
        }
        if ( ! has_filter ) { continue; }
        // check filter logic
        if ( filt != NULL ) {
            int32_t ret = filter_test(filt, iv, NULL);
            if ( filter_logic == FLT_INCLUDE ) { if ( !ret)  has_filter = false; }
            else if ( ret ) { has_filter = false; }
        }
        if ( ! has_filter ) { continue; }

        if ( bcf_get_format_int32(odr.hdr, iv, "DP", &v_dp, &n_fmt) < 0 ) {
            error("[E:%s:%d %s] Cannot find the field DP from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
        }
        n_fmt = 0;
        if ( bcf_get_format_int32(odr.hdr, iv, "RBQ", &v_rbq, &n_fmt) < 0 ) {continue;}
        n_fmt = 0;
        if ( bcf_get_format_int32(odr.hdr, iv, "ABQ", &v_abq, &n_fmt) < 0 ) {continue;}
        n_fmt = 0;
        if ( bcf_get_format_int32(odr.hdr, iv, "RDF", &v_rdf, &n_fmt) < 0 ) {continue;}
        n_fmt = 0;
        if ( bcf_get_format_int32(odr.hdr, iv, "ADF", &v_adf, &n_fmt) < 0 ) {continue;}
        n_fmt = 0;
        if ( bcf_get_format_int32(odr.hdr, iv, "RDR", &v_rdr, &n_fmt) < 0 ) {continue;}
        n_fmt = 0;
        if ( bcf_get_format_int32(odr.hdr, iv, "ADR", &v_adr, &n_fmt) < 0 ) {continue;}
        n_fmt = 0;

        // Decide which strand to use (?)
        if (single_strand) {
            int32_t DPF = 0, DPR = 0;
            for (int32_t i=0; i < nsamples; ++i) {
                if ( v_dp[i] == bcf_int32_missing ) { continue; }
                DPF += v_rdf[i] + v_adf[i];
                DPR += v_rdr[i] + v_adr[i];
            }
            int32_t strand = (DPF > DPR) ? 1 : 0;
            for (int32_t i=0; i < nsamples; ++i) {
                if ( v_dp[i] == bcf_int32_missing ) { continue; }
                if ( strand == 0 ) {
                    v_rdf[i] = v_rdr[i];
                    v_adf[i] = v_adr[i];
                }
            }
        } else {
            for (int32_t i=0; i < nsamples; ++i) {
                v_rdf[i] = v_rdf[i] + v_rdr[i];
                v_adf[i] = v_adf[i] + v_adr[i];
            }
        }

        int32_t max_ctrl_rdp = 0;
        int32_t max_ctrl_adp = 0;
        double  max_raf_ctrl = 0;
        double  max_aaf_ctrl = 0;
        int32_t obs_ctrl = 0;
        int32_t obs_test = 0;
        int32_t max_sum_ctrl_rbq = 0;
        int32_t max_sum_ctrl_abq = 0;
        int32_t max_sum_test_rbq = 0;
        int32_t max_sum_test_abq = 0;
        double  sum_test_rdp = 0;
        double  sum_test_adp = 0;
        int32_t ctrl_homo = -1;
        int32_t sum_ctrl_rdp = 0;
        int32_t sum_ctrl_adp = 0;

        std::string vtype = bcf_is_snp(iv) ? "SNP" : "INDEL";
        std::stringstream ss_indiv;
        for (auto & i : test_sample_idx) {
            if ( v_dp[i] == bcf_int32_missing || v_dp[i] < minDP ) {
                if (obs_test > 0) {ss_indiv << "\t";}
                ss_indiv << "-1,-1";
                continue;
            }
            sum_test_rdp += v_rdf[i];
            sum_test_adp += v_adf[i];
            if (max_sum_test_rbq < v_rbq[i] * v_rdf[i]) {
                max_sum_test_rbq = v_rbq[i] * v_rdf[i];
            }
            if (max_sum_test_abq < v_abq[i] * v_adf[i]) {
                max_sum_test_abq = v_abq[i] * v_adf[i];
            }
            if (obs_test > 0) {ss_indiv << "\t";}
            ss_indiv << v_rdf[i] << "," << v_adf[i];
            obs_test++;
        }
        for(auto & i : ctrl_sample_idx) {
            if ( v_dp[i] == bcf_int32_missing || v_dp[i] < minDP ) {
                ss_indiv << "\t-1,-1";
                continue;
            }
            obs_ctrl++;
            sum_ctrl_rdp += v_rdf[i];
            sum_ctrl_adp += v_adf[i];
            max_ctrl_rdp = max_ctrl_rdp < v_rdf[i] ? v_rdf[i] : max_ctrl_rdp;
            max_ctrl_adp = max_ctrl_adp < v_adf[i] ? v_adf[i] : max_ctrl_adp;
            int32_t rbq_sum = v_rbq[i] * v_rdf[i];
            int32_t abq_sum = v_abq[i] * v_adf[i];
            double  vaf = (double)v_adf[i] / (double)(v_rdf[i] + v_adf[i]);
            if (max_sum_ctrl_rbq < rbq_sum) {
                max_sum_ctrl_rbq = rbq_sum;
            }
            if (max_sum_ctrl_abq < abq_sum) {
                max_sum_ctrl_abq = abq_sum;
            }
            if (max_raf_ctrl < 1.-vaf) {
                max_raf_ctrl = 1.-vaf;
            }
            if (max_aaf_ctrl < vaf) {
                max_aaf_ctrl = vaf;
            }
            ss_indiv << "\t" << v_rdf[i] << "," << v_adf[i];
        }
        if (obs_ctrl < 1 || obs_test < 1) {
            continue;
        }
        if (obs_ctrl >= min_ctrl_support) {
            if (max_raf_ctrl < max_minor_vaf_ctrl &&
                max_ctrl_rdp < max_minor_adp_ctrl &&
                max_sum_ctrl_rbq < max_sum_minor_abq_ctrl) {
                ctrl_homo = 1;
            } else if (max_aaf_ctrl < max_minor_vaf_ctrl &&
                max_ctrl_adp < max_minor_adp_ctrl &&
                max_sum_ctrl_abq < max_sum_minor_abq_ctrl) {
                ctrl_homo = 0;
            }
        }
        std::stringstream ss;
        ss << ctrl_homo << "\t";
        ss << obs_ctrl << ":" << max_ctrl_rdp << "," << max_ctrl_adp << "\t";
        ss << obs_test << ":" << max_sum_test_rbq << "," << max_sum_test_abq;
        ++nVariant;

        if (preview && ctrl_homo != -1) {
            double test_vaf = (double) sum_test_adp / (sum_test_adp + sum_test_rdp);
            if ((ctrl_homo == 0 && sum_test_adp >= min_minor_adp_test && test_vaf >= min_minor_vaf_test) ||
                (ctrl_homo == 1 && sum_test_rdp >= min_minor_adp_test && 1 - test_vaf > min_minor_vaf_test)) {
                printf("%s\t%ld\t%s\t%s\t%s\t%s\n", bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, vtype.c_str(), bcf_hdr_int2id(odr.hdr,BCF_DT_ID,iv->d.flt[0]), ss.str().c_str(), ss_indiv.str().c_str() );
            }
        }

        hprintf(wf, "%s\t%d\t%s\t%s\t%s\n", bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, vtype.c_str(), ss.str().c_str(), ss_indiv.str().c_str() );
    }
    hts_close(wf);
    return 0;
}
