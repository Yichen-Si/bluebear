#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "hts_utils.h"
#include <bits/stdc++.h>
#include <sstream>
#include "utils.h"

int32_t cmdVcfCarrier(int32_t argc, char** argv) {

	std::string inVcf, reg, out;
	int32_t min_ac=1, max_ac=50;
	int32_t verbose = 10000;
	bool snp_only = false;
  	bool haploid = false;

	bcf_vfilter_arg vfilt;

	paramList pl;
	BEGIN_LONG_PARAMS(longParameters)
		LONG_PARAM_GROUP("Input", NULL)
		LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
		LONG_STRING_PARAM("region",&reg, "Region")

		LONG_PARAM_GROUP("Variant Filtering Options", NULL)
		LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
		LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
		LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")
		LONG_PARAM("snp-only",&snp_only, "")
		LONG_PARAM("haploid",&haploid,"Input is haploid though vcf is read as diploid (?)")

		LONG_PARAM_GROUP("Additional Options", NULL)
		LONG_INT_PARAM("min-ac",&min_ac, "Minimum sample allele count to consider")
		LONG_INT_PARAM("max-ac",&max_ac, "Maximum sample allele count to consider")
		LONG_INT_PARAM("verbose",&verbose, "Periodic report when scanning VCF")

		LONG_PARAM_GROUP("Output Options", NULL)
		LONG_STRING_PARAM("out", &out, "Output file prefix")

	END_LONG_PARAMS();
	pl.Add(new longParams("Available Options", longParameters));
	pl.Read(argc, argv);
	pl.Status();

	if (inVcf.empty() || out.empty()) {
		error("[E:%s:%d %s] --in-vcf and --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
	}

	// bcf reader
	std::vector<GenomeInterval> intervals;
	if (!reg.empty())
		parse_intervals(intervals, "", reg);
	BCFOrderedReader odr(inVcf, intervals);
	bcf1_t* iv = bcf_init();
	int32_t n_samples = bcf_hdr_nsamples(odr.hdr);
	char**  id_ptr = bcf_hdr_get_samples(odr.hdr);
	std::vector<std::string> sample_id(id_ptr, id_ptr+n_samples);
	std::map<std::string, int32_t> id_map; // sample_id -> index in genotype vector
	for (auto i = 0; i < n_samples; ++i) {
		id_map.insert({sample_id[i], i});
	}
    int32_t nalleles = 2*n_samples;
    if (haploid) {nalleles = n_samples;}
	if ( n_samples == 0 )
		error("FATAL ERROR: The VCF does not have any samples with genotypes");
	notice("VCF contain %d samples", n_samples);

	// output
	std::string outf = out + ".ac_" + std::to_string(min_ac) + "_" + std::to_string(max_ac) + ".tsv.gz";
	htsFile* wf = hts_open(outf.c_str(), "wz");
	hprintf(wf, "#CHROM\tPOS\tALLELE\tAC\tMAJOR\tMINOR\tFILTER\tCARRIERS\n");

	// handle filter string
	vfilt.init(odr.hdr);

	int32_t* p_gt = NULL;
	int32_t* info_ac = NULL;
	int32_t* info_an = NULL;
	int32_t n_gt = 0, n_ac = 0, n_an = 0;
	int32_t nVariant = 0;

	std::string major, minor;

	for(int32_t k=0; odr.read(iv); ++k) {

		if (!bcf_is_snp(iv) && snp_only) {continue;}

		bcf_unpack(iv, BCF_UN_STR);
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

		// Allele info
		if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
		if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
        int32_t ac = info_ac[0];
		bool flip = 0;
		std::string ref(iv->d.allele[0]);
		std::string alt(iv->d.allele[1]);
        if (info_an[0] > nalleles && haploid) {ac = ac / 2;}
        if (ac > nalleles / 2) {
			flip = 1;
			major = alt;
			minor = ref;
			ac = nalleles - ac;
		} else {
			major = ref;
			minor = alt;
		}
        if (ac < min_ac || ac > max_ac) {continue;}

		// extract genotype
		if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
			error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
		}

		// output
		// CHROM, POS, ALLELE, AC, MAJOR, MINOR,
		hprintf(wf, "%s\t%ld\t%s,%s\t%d\t%s\t%s\t", bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, ref.c_str(), alt.c_str(), ac, major.c_str(), minor.c_str());
		// FILTER
		if (iv->d.n_flt == 0) {
			hprintf(wf, ".\t");
		}
		else {
			hprintf(wf, "%s", bcf_hdr_int2id(odr.hdr, BCF_DT_ID, iv->d.flt[0]) );
			if (iv->d.n_flt > 1) {
				for(int32_t i=1; i < iv->d.n_flt; ++i) {
					hprintf(wf, ",%s", bcf_hdr_int2id(odr.hdr, BCF_DT_ID, iv->d.flt[i]) );
				}
			}
			hprintf(wf, "\t");
		}
		// CARRIERS
		for (int32_t i = 0; i < n_samples; ++i) {
			if ( bcf_gt_is_missing(p_gt[2*i]) || bcf_gt_is_missing(p_gt[2*i+1]) ) {
				continue;
			}
			int32_t g1 = bcf_gt_allele(p_gt[2*i]);
			int32_t g2 = bcf_gt_allele(p_gt[2*i+1]);
			int32_t geno;
			if (flip) {
				if (g1 > 0 && g2 > 0) {continue;}
				g1 = 1 - g1;
				g2 = 1 - g2;
			}
			if (g1 + g2 >= 1) { // Contain minor allele
				hprintf(wf, "%s,%d%d;", sample_id[i].c_str(), g1, g2);
			}
		}
		hprintf(wf, "\n");
		nVariant ++;
		// periodic message to user
		if ( k % verbose == 0 )
			notice("Processing %d markers at %s:%d. Recorded %d rare variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nVariant);
	}
	hts_close(wf);
	return 0;
}
