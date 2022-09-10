#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include <bits/stdc++.h>
#include <sstream>

#include "utils.h"
#include "pbwt_build.h"
#include "parsimony_tree.h"

int32_t temp(int32_t argc, char** argv) {

	// std::string inVcf, reg, out;
	// int32_t verbose = 10000;
	// bool haploid = 1;

	// bcf_vfilter_arg vfilt;
	// bcf_gfilter_arg gfilt;

	// paramList pl;
	// BEGIN_LONG_PARAMS(longParameters)
	// 	LONG_PARAM_GROUP("Input", NULL)
	// 	LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
	// 	LONG_STRING_PARAM("region",&reg, "Region")

	// 	LONG_PARAM_GROUP("Variant Filtering Options", NULL)
	// 	LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
	// 	LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
	// 	LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

	// 	LONG_PARAM_GROUP("Additional Options", NULL)
	// 	LONG_INT_PARAM("verbose",&verbose, "Periodic report when scanning VCF")
	// 	LONG_INT_PARAM("haploid",&haploid, "0/1 for whether to treat input as haploid")

	// 	LONG_PARAM_GROUP("Output Options", NULL)
	// 	LONG_STRING_PARAM("out", &out, "Output file prefix")

	// END_LONG_PARAMS();
	// pl.Add(new longParams("Available Options", longParameters));
	// pl.Read(argc, argv);
	// pl.Status();

	// if (inVcf.empty() || out.empty()) {
	// 	error("[E:%s:%d %s] --in-vcf and --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
	// }

	// // bcf reader
	// std::vector<GenomeInterval> intervals;
	// if (!reg.empty())
	// 	parse_intervals(intervals, "", reg);
	// BCFOrderedReader odr(inVcf, intervals);
	// bcf1_t* iv = bcf_init();
	// int32_t n_samples = bcf_hdr_nsamples(odr.hdr);
	// if ( n_samples == 0 )
	// 	error("FATAL ERROR: The VCF does not have any samples with genotypes");
	// notice("VCF contain %d samples", n_samples);

	// // output
	// // std::string outf = out + ".";
	// // std::ofstream wf(outf.c_str(), std::ios::out);
	// // wf << "#AC\tPOS\tALLELE\tMAJOR\tMINOR\tFILTER\tCARRIERS\n";

	// // handle filter string
	// std::string filter_str;
	// int32_t filter_logic = 0;
	// if ( vfilt.include_expr.empty() ) {
	// 	if ( vfilt.exclude_expr.empty() ) {
	// 	} else {
	// 		filter_str = vfilt.exclude_expr;
	// 		filter_logic |= FLT_EXCLUDE;
	// 	}
	// } else {
	// 	if ( vfilt.exclude_expr.empty() ) {
	// 		filter_str = vfilt.include_expr;
	// 		filter_logic |= FLT_INCLUDE;
	// 	} else {
	// 		error("[E:%s:%d %s] Cannot use both --include-expr and --exclude-expr options",__FILE__,__LINE__,__FUNCTION__);
	// 	}
	// }

	// filter_t* filt = NULL;
	// if ( filter_logic != 0 )
	// 	filter_init(odr.hdr, filter_str.c_str());

	// // handle --apply-filtrs
	// std::vector<int32_t> req_flt_ids;
	// if ( !vfilt.required_filters.empty() ) {
	// 	for(int32_t i=0; i < (int32_t)vfilt.required_filters.size(); ++i) {
	// 		req_flt_ids.push_back(bcf_hdr_id2int(odr.hdr, BCF_DT_ID, vfilt.required_filters[i].c_str()));
	// 	}
	// }

	// parsimonyTree pTree;

	// int32_t *p_gt, *info_ac, *info_an;
	// p_gt = info_ac = info_an = NULL;
	// int32_t n_gt = 0, n_ac = 0, n_an = 0;
	// int32_t nVariant = 0;

	// for(int32_t k=0; odr.read(iv); ++k) {

	// 	if (!bcf_is_snp(iv) && snp_only) {continue;}
	// 	// unpack FILTER column
	// 	bcf_unpack(iv, BCF_UN_STR);
	// 	// check --apply-filters
	// 	bool has_filter = req_flt_ids.empty() ? true : false;
	// 	if ( ! has_filter ) {
	// 		for(int32_t i=0; i < iv->d.n_flt; ++i) {
	// 			for(int32_t j=0; j < (int32_t)req_flt_ids.size(); ++j) {
	// 				if ( req_flt_ids[j] == iv->d.flt[i] )
	// 					has_filter = true;
	// 			}
	// 		}
	// 	}
	// 	if ( ! has_filter ) { continue; }
	// 	// check filter logic
	// 	if ( filt != NULL ) {
	// 		int32_t ret = filter_test(filt, iv, NULL);
	// 		if ( filter_logic == FLT_INCLUDE ) { if ( !ret)  has_filter = false; }
	// 		else if ( ret ) { has_filter = false; }
	// 	}
	// 	if ( ! has_filter ) { continue; }

	// 	// periodic message to user
	// 	if ( k % verbose == 0 )
	// 		notice("Processed %d variants, current at %s:%d.", nVariant, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);

	// 	// Allele info
	// 	if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
	// 	if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
	// 	std::string ref(iv->d.allele[0]);
	// 	std::string alt(iv->d.allele[1]);
	// 	// extract genotype
	// 	if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
	// 		error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
	// 	}
	// 	// TODO: now assume haploid and use only the first allele
	// 	// TODO: need to handle missingness
	// 	if (pTree.add_variable_site_vcf(iv) ) {
	// 		nVariant ++;
	// 	}
	// 	// bool gt[n_samples];
	// 	// for (int32_t i = 0; i < n_samples; ++i) {
	// 	// 	int32_t g1 = p_gt[2*i];
	// 		// int32_t g2 = p_gt[2*i+1];
	// 		// int32_t geno;
	// 		// if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
	// 		// 	continue;
	// 		// } else {
	// 		// 	geno = (int32_t) ((bcf_gt_allele(g1) > 0) && (bcf_gt_allele(g2) > 0));
	// 		// }
	// 	// 	gt[i] = bcf_gt_is_missing(g1);
	// 	// 	if (! gt[i])  {
	// 	// 		gt[i] = bcf_gt_allele(g1) > 0;
	// 	// 	}
	// 	// }

	// }

	// if (pTree.build_pbwt() > 0) {
	// 	error("PBWT failed\n");
	// }
	// if (pTree.run_greedy_parsimony()) {
	// 	pTree.clearn_rebuild()
	// }





	return 0;
}
