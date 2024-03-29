#include "utils.h"
#include "cramore.h"
#include "hts_utils.h"

#include <fstream>

#include "compact_matrix.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"


int32_t temp(int32_t argc, char** argv) {

  std::string inVcf, reg, out, samples, info, vmap;
  int32_t verbose = 10000;
  double minDprime = -1;
  bool haploid = false;
  bool count_only = false;
  bool snp_only = false;
  bool anno_group = false;
  bool debug = false;

  bcf_vfilter_arg vfilt;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
	LONG_PARAM_GROUP("Input Sites", NULL)
	LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
	LONG_STRING_PARAM("samples",&samples,"List of sample ID")
    LONG_STRING_PARAM("sample-info",&info,"Sample population or other group info")
	LONG_STRING_PARAM("variant-map",&vmap,"Mapping between individual variant and collapsed ID")

	LONG_PARAM_GROUP("Variant Filtering Options", NULL)
	LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
	LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
	LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")
    LONG_PARAM("snp-only",&snp_only,"Exclude indels")
    LONG_PARAM("haploid",&haploid,"Input is haploid though vcf is read as diploid (?)")

	LONG_PARAM_GROUP("Output Options", NULL)
	LONG_STRING_PARAM("out", &out, "Output file prefix")
	LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")
    LONG_PARAM("debug",&debug,"")



  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

    // sanity check of input arguments
    if ( inVcf.empty() || out.empty() || vmap.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
    }

    struct VarConfig {
        int32_t ac, n;
        std::map<std::string, int32_t> group_ct;
        std::map<int32_t, int32_t> filter_ct;
        std::string first_occur;
        std::string id;
    };

    std::vector<GenomeInterval> intervals;
    if ( !reg.empty() ) {
        parse_intervals(intervals, "", reg);
    }
    BCFOrderedReader odr(inVcf, intervals);
    bcf1_t* iv = bcf_init();
    int32_t pass_id = (int32_t) bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "PASS");
    if (!samples.empty()) {
        bcf_hdr_set_samples(odr.hdr, samples.c_str(), 1);
    }

    // Output collapsed variant configurations
    std::string outf = out + ".var.config.tsv.gz";
    htsFile* wf = hts_open(outf.c_str(), "wg");
    std::map<int32_t, std::map<std::string, VarConfig> > clpsVars;

	// Read variant mapping
	std::fstream vf;
	vf.open(vmap.c_str(), std::fstream::in);
	if (!vf.is_open()) {
		error("Unable to read provided variant map info");
	}
	std::string line;
	std::vector<std::string> words;
	std::map<std::string, std::pair<int32_t, std::string> > config_map;
	int32_t n_missing = 0;
	while(std::getline(vf, line)) {
		split(words, "\t", line);
		if (words.size() < 2) {
			continue;
		}
		config_map[words[0]] = std::pair<int32_t, std::string> {std::stoi(words[1]), words[2]};
	}

    // Read sample ID
    int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
	char** id_ptr = bcf_hdr_get_samples(odr.hdr);
	std::vector<std::string> sample_id(id_ptr, id_ptr+nsamples);
    std::map<int32_t, std::string> sample_info;
    std::map<std::string, int32_t> sample_indx;
    for (uint32_t i = 0; i < nsamples; ++i) {
        sample_indx[sample_id[i]] = i;
        sample_info[i] = "Unknown";
    }
    notice("Detected %d samples", nsamples);
    // Read sample group info
    std::set<std::string> group_set;
    std::vector<std::string> group_list;
    if (!info.empty()) {
        std::fstream rf;
        rf.open(info.c_str(), std::fstream::in);
        if (!rf.is_open()) {
            error("Unable to read provided sample group info");
        }
        words.clear();
        int32_t n_missing = 0;
        while(std::getline(rf, line)) {
            split(words, "\t", line);
            if (words.size() < 2) {
                continue;
            }
            group_set.insert(words[1]);
            auto ptr = sample_indx.find(words[0]);
            if (ptr != sample_indx.end()) {
                sample_info[ptr->second] = words[1];
            }
        }
        for (uint32_t i = 0; i < nsamples; ++i) {
            if (sample_info[i] == "Unknown") {
                n_missing += 1;
            }
        }
        group_list.resize(group_set.size());
        std::copy(group_set.begin(), group_set.end(), group_list.begin());
        std::sort(group_list.begin(), group_list.end());
        if (n_missing > 0) {
            group_list.push_back("Unknown");
        }
        notice("Missing group info for %d individuals", n_missing);

        if (n_missing < nsamples)
            anno_group = true;
    }

    // Output annotated variant info
    outf = out + ".site.vcf.gz";
    BCFOrderedWriter odw(outf.c_str(),0);
    bcf_hdr_t* hnull = bcf_hdr_subset(odr.hdr, 0, 0, 0);
    bcf_hdr_remove(hnull, BCF_HL_FMT, NULL);
    odw.set_hdr(hnull);
	// Add info field
	char buffer[65536];
	// check the existence of header and create one if needed
	if ( anno_group && bcf_hdr_id2int(odw.hdr, BCF_DT_ID, "GroupAC") < 0 ) {
		sprintf(buffer,"##INFO=<ID=GroupAC,Number=1,Type=String,Description=\"Minor allele count by group/population\">\n");
		bcf_hdr_append(odw.hdr, buffer);
	}
	if ( bcf_hdr_id2int(odw.hdr, BCF_DT_ID, "ConfigID") < 0 ) {
		sprintf(buffer,"##INFO=<ID=ConfigID,Number=1,Type=String,Description=\"ID of the collapsed unique variant\">\n");
		bcf_hdr_append(odw.hdr, buffer);
	}
    odw.write_hdr();

    if (!samples.empty()) {
        bcf_hdr_set_samples(odr.hdr, samples.c_str(), 1);
    }

    // handle filter string
    vfilt.init(odr.hdr);
    vfilt.snpOnly = snp_only;
    int32_t nalleles = 2*nsamples;
    if (haploid) {nalleles = nsamples;}

    int32_t  n_gt = 0, n_ac = 0, n_an = 0;
    int32_t* p_gt = NULL;
    int32_t  nVariant = 0, nUniq = 0;
    int32_t* info_ac = NULL;
    int32_t* info_an = NULL;

    for(int32_t k=0; odr.read(iv); ++k) {  // read marker

        if ( k % verbose == 0 )
            notice("Processing %d markers at %s:%d, recorded %d unique configurations.", nVariant, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nUniq);

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

        std::stringstream ss;
		ss << iv->pos+1 << "_" <<  bcf_get_allele(iv)[0] << "_" << bcf_get_allele(iv)[1];
		std::string variant_id = ss.str();
		auto v_ptr = config_map.find(variant_id);
		if (v_ptr == config_map.end()) {continue;}

		int32_t ac = (v_ptr->second).first;
		std::string config_id = (v_ptr->second).second;

        std::map<std::string, int32_t> group_count;
        for (const auto & v : group_list) {
            group_count[v] = 0;
        }
        int32_t ct = 0;
        auto it = clpsVars.find(ac);
        if (it == clpsVars.end()) {
            clpsVars[ac] = std::map<std::string, VarConfig> {};
        }
        it = clpsVars.find(ac);
        auto ptr = (it->second).find(config_id);
		if (ptr == (it->second).end()) {
			if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
				error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
			}
            if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
            if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
            bool flip = info_ac[0] > info_an[0] / 2;
			if (haploid) {
				for(int32_t i=0; i < nsamples; ++i) {
					uint8_t geno = 0;
					if ( bcf_gt_is_missing(p_gt[i*2]) ) {
						geno = 0;
					} else {
						geno = (bcf_gt_allele(p_gt[i*2]) > 0) ? 1 : 0;
					}
					if (flip) {
						geno = 1 - geno;
					}
					group_count[sample_info[i]] += geno;
				}
			} else {
				for(int32_t i=0; i < nsamples; ++i) {
					uint8_t geno = 0;
					if ( bcf_gt_is_missing(p_gt[i*2]) | bcf_gt_is_missing(p_gt[i*2+1]) ) {
						geno = 0;
					} else {
						geno = ((bcf_gt_allele(p_gt[i*2]) > 0) ? 1 : 0) + ((bcf_gt_allele(p_gt[i*2+1]) > 0) ? 1 : 0);
					}
					if (flip) {
						geno = 2 - geno;
					}
					group_count[sample_info[i]] += geno;
				}
			}
            VarConfig new_var;
            new_var.ac = ac;
            new_var.n  = 1;
            if (anno_group)
                new_var.group_ct = group_count;
            for(int32_t i=0; i < iv->d.n_flt; ++i)
                new_var.filter_ct[iv->d.flt[i]] = 1;
            new_var.id = config_id;
            new_var.first_occur = variant_id;
            (it->second)[config_id] = new_var;
            nUniq += 1;
        } else {
            (ptr->second).n += 1;
            for(int32_t i=0; i < iv->d.n_flt; ++i)
                (ptr->second).filter_ct[iv->d.flt[i]] += 1;
        }

        ++nVariant;

        bcf1_t* nv = bcf_dup(iv);
        bcf_unpack(nv, BCF_UN_ALL);
        bcf_subset(odw.hdr, nv, 0, 0);

        if (anno_group) {
            ss.str(std::string());
            for (const auto & v : group_list) {
                ss << v << ":" << group_count[v] << ".";
            }
            bcf_update_info_string(odw.hdr, nv, "GroupAC", ss.str().c_str());
        }
        bcf_update_info_string(odw.hdr, nv, "ConfigID", config_id.c_str());
        odw.write(nv);
        bcf_destroy(nv);
        if (debug && nUniq > 50) {
            break;
        }
    }
    odr.close();

    notice("Finish processing all variants, identified %d unique ones.\n", nUniq);

    for (auto & w : clpsVars) {
        for (auto & v : w.second) {
            // ID AC nVar nPASS nFAIL FAILdetail CarrierPOP Config
            int32_t nPASS = v.second.filter_ct[pass_id];
            int32_t nFAIL = v.second.n - nPASS;
            std::stringstream ss;
            for (const auto & w : v.second.filter_ct) {
                ss << bcf_hdr_int2id(odw.hdr,BCF_DT_ID,w.first) << ":" << w.second << ".";
            }
            std::string filter_detail = ss.str();
            ss.str(std::string());
            for (const auto & w : v.second.group_ct) {
                ss << w.first << ":" << w.second << ".";
            }
            std::string carrier_group = ss.str();

            hprintf(wf, "%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n", v.second.id.c_str(), v.second.ac, v.second.n, nPASS, nFAIL, filter_detail.c_str(), carrier_group.c_str(), v.second.first_occur.c_str());
        }
    }
    odw.close();
    hts_close(wf);
    return 0;
}
