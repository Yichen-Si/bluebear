#include "cramore.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "hts_utils.h"
#include <bits/stdc++.h>
#include <sstream>
#include "utils.h"

int32_t cmdVcfAnnoHapGroup(int32_t argc, char** argv) {

	std::string inVcf, inGroup, out, reg;
	int32_t min_ac=2, mlevel=0, olevel=0;
	int32_t verbose = 10000;

	bcf_vfilter_arg vfilt;
	bcf_gfilter_arg gfilt;

	paramList pl;
	BEGIN_LONG_PARAMS(longParameters)
		LONG_PARAM_GROUP("Input", NULL)
		LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
		LONG_STRING_PARAM("region",&reg, "Region")
		LONG_STRING_PARAM("in-group",&inGroup,"Input sample and group info")

		LONG_PARAM_GROUP("Variant Filtering Options", NULL)
		LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
		LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
		LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")

		LONG_PARAM_GROUP("Additional Options", NULL)
		LONG_INT_PARAM("min-ac",&min_ac, "Minimum sample allele count to consider")
		LONG_INT_PARAM("verbose",&verbose, "Periodic report when scanning VCF")
		LONG_INT_PARAM("max-level",&mlevel, "Maximum level of the hierarchical group to use")

		LONG_PARAM_GROUP("Output Options", NULL)
		LONG_STRING_PARAM("out", &out, "Output file prefix")
		LONG_INT_PARAM("out-level",&olevel, "")

	END_LONG_PARAMS();
	pl.Add(new longParams("Available Options", longParameters));
	pl.Read(argc, argv);
	pl.Status();

	if (inVcf.empty() || inGroup.empty() || out.empty()) {
		error("[E:%s:%d %s] --in-vcf, --in-group, and --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
	}

	// bcf reader
	std::vector<GenomeInterval> intervals;
	if (!reg.empty())
		parse_intervals(intervals, "", reg);
	BCFOrderedReader odr(inVcf, intervals);
	bcf1_t* iv = bcf_init();
	int32_t n_samples = bcf_hdr_nsamples(odr.hdr);
	char** id_ptr = bcf_hdr_get_samples(odr.hdr);
	std::vector<std::string> sample_id(id_ptr, id_ptr+n_samples);
	std::map<std::string, int32_t> id_map; // sample_id -> index in genotype vector
	for (auto i = 0; i < n_samples; ++i) {
		id_map.insert({sample_id[i], i});
	}
	if ( n_samples == 0 )
		error("FATAL ERROR: The VCF does not have any samples with genotypes");
	notice("VCF contain %d samples", n_samples);

	// Read sample group file
	std::ifstream rf;
	std::string line;
	std::vector<std::string> v;
	std::set<std::string> branch_set;
	std::map<std::string, std::vector<std::string> > branch_level; // full code -> [level0, level1, ...]
	std::map<std::string, std::string> id_info;
	rf.open(inGroup, std::ifstream::in);
	while (std::getline(rf, line)) {
			split(v, "\t", line);
			if ((int32_t) v.size() < mlevel+2) {
					error("Invalid --in-group.");
			}
			if (id_map.find(v[0]) == id_map.end()) {
				continue;
			}
			id_info.insert({v[0], v[1]});
			if ( branch_level.find(v[1]) == branch_level.end() ) {
				std::vector<std::string> tv (v.begin()+1, v.end());
				branch_level.insert({v[1], tv});
			}
			branch_set.insert(v[1]);
	}
	rf.close();

// for (auto & x : branch_level) {
//   std::cout << x.first;
//   for (auto & y : x.second) {
//     std::cout << '\t' << y;
//   }
//   std::cout << '\n';
// }

	int32_t nbranch = branch_set.size();
	if (olevel > mlevel) {
		olevel = mlevel;
	}
	std::vector<std::string> branch(branch_set.begin(), branch_set.end());
	std::sort(branch.begin(), branch.end());
	auto itr = id_map.begin();
	while (itr != id_map.end()) {
			if (id_info.find(itr->first) == id_info.end()) {
				 itr = id_map.erase(itr);
			} else {
				 ++itr;
			}
	}
	std::map<std::string, int32_t> branch_map;
	for (int32_t i = 0; i < nbranch; ++i) {
		branch_map.insert({branch[i], i});
	}
	notice("Read group info, Total %d branches", nbranch);

	// output
	std::string out_table = out + ".count.cross_only.all_branch.tsv";
	std::ofstream wf(out_table.c_str(), std::ios::out);
	wf << "POS\tREF\tALT\tQTAL\tFILTER\tINFO";
	for (auto & v : branch) {
		wf << '\t' << v;
	}
	std::string out_vcf = out + ".sites.vcf.gz";
	BCFOrderedWriter odw(out_vcf.c_str(),0);
	bcf_hdr_t* hnull = bcf_hdr_subset(odr.hdr, 0, 0, 0);
	bcf_hdr_remove(hnull, BCF_HL_FMT, NULL);
	odw.set_hdr(hnull);

	// Add info field
	char buffer[65536];
	// check the existence of header and create one if needed
	if ( bcf_hdr_id2int(odw.hdr, BCF_DT_ID, "n_group_REF") < 0 ) {
		sprintf(buffer,"##INFO=<ID=n_group_REF,Number=.,Type=Integer,Description=\"\">\n");
		bcf_hdr_append(odw.hdr, buffer);
	}
	if ( bcf_hdr_id2int(odw.hdr, BCF_DT_ID, "n_group_ALT") < 0 ) {
		sprintf(buffer,"##INFO=<ID=n_group_ALT,Number=.,Type=Integer,Description=\"\">\n");
		bcf_hdr_append(odw.hdr, buffer);
	}
	if ( bcf_hdr_id2int(odw.hdr, BCF_DT_ID, "YCCminor_1st") < 0 ) {
		sprintf(buffer,"##INFO=<ID=YCCminor_1st,Number=.,Type=String,Description=\"\">\n");
		bcf_hdr_append(odw.hdr, buffer);
	}
	if ( bcf_hdr_id2int(odw.hdr, BCF_DT_ID, "YCCminor_2nd") < 0 ) {
		sprintf(buffer,"##INFO=<ID=YCCminor_2nd,Number=.,Type=String,Description=\"\">\n");
		bcf_hdr_append(odw.hdr, buffer);
	}
	odw.write_hdr();

	// handle filter string
	std::string filter_str;
	int32_t filter_logic = 0;
	if ( vfilt.include_expr.empty() ) {
		if ( vfilt.exclude_expr.empty() ) {
		} else {
			filter_str = vfilt.exclude_expr;
			filter_logic |= FLT_EXCLUDE;
		}
	} else {
		if ( vfilt.exclude_expr.empty() ) {
			filter_str = vfilt.include_expr;
			filter_logic |= FLT_INCLUDE;
		} else {
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

	int32_t* p_gt = NULL;
	int32_t* info_ac = NULL;
	int32_t* info_an = NULL;
	int32_t n_gt = 0, n_ac = 0, n_an = 0;
	int32_t nVariant = 0;

	for(int32_t k=0; odr.read(iv); ++k) {
		// periodic message to user
		if ( k % verbose == 0 )
			notice("Processing %d markers at %s:%d. Recorded %d rare variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nVariant);
		// unpack FILTER column
		bcf_unpack(iv, BCF_UN_STR);
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

		// extract genotype and apply genotype level filter
		if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
			error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
		}
		if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
		if (info_ac[0] < min_ac) {continue;}
		if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
		if (info_an[0]-info_ac[0] < min_ac) {continue;}
		bcf1_t* nv = bcf_dup(iv);
		bcf_subset(odw.hdr, nv, 0, 0);
		bcf_unpack(nv, BCF_UN_STR);

		std::vector<int32_t> count_ref (nbranch, 0);
		std::vector<int32_t> count_alt (nbranch, 0);
		std::vector<std::string> topycc{"", ""};

		for (const auto & v : id_map) {
			int32_t i = v.second;
			int32_t g1 = p_gt[2*i];
			int32_t g2 = p_gt[2*i+1];
			int32_t geno;
			if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
				continue;
			} else {
				geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
				if (geno >= 1) { // ALT
					count_alt[branch_map[id_info[v.first]]] += 1;
				} else {
					count_ref[branch_map[id_info[v.first]]] += 1;
				}
			}
		} // Count alleles by group
		// Number of different group for each allele at each level
		std::vector<int32_t> ngr(mlevel+1, 0);
		std::vector<int32_t> nga(mlevel+1, 0);
		for (int32_t i = 0; i < nbranch; ++i) {
			if (count_ref[i] > 0) {
				ngr[0]++;
			}
			if (count_alt[i] > 0) {
				nga[0]++;
			}
		}

		std::vector<std::vector<std::pair<std::string, int32_t> > > top_groups;
		int32_t derived = 1; // Guess which allele is the derived allele
		if (ngr[0] < nga[0]) {
			derived = 0; // Reference allele is more restricted on the tree
		}
		std::vector<std::pair<std::string, int32_t> > top_pair;
		std::vector<int32_t> top_index(count_ref.size());
		std::iota(top_index.begin(), top_index.end(), 0);
		if (derived == 0) {
			sort(top_index.begin(), top_index.end(),
				[&count_ref](int32_t left, int32_t right) -> bool {
						return count_ref[left] > count_ref[right]; } );
			for (int32_t rr = 0; rr < 2; ++rr) {
				std::stringstream sstr;
				int32_t r = top_index[rr];
				top_pair.push_back({branch[r], count_ref[r] } );
				if (count_ref[r] == 0) {
					sstr << "-,0.";
				} else {
					sstr << branch[r] << "," << count_ref[r] << ".";
				}
				topycc[rr] += sstr.str();
			}
		} else {
			sort(top_index.begin(), top_index.end(),
				[&count_alt](int32_t left, int32_t right) -> bool {
						return count_alt[left] > count_alt[right]; } );
			for (int32_t rr = 0; rr < 2; ++rr) {
				int32_t r = top_index[rr];
				top_pair.push_back({branch[r], count_alt[r] } );
				std::stringstream sstr;
				if (count_alt[r] == 0) {
					sstr << "-,0.";
				} else {
					sstr << branch[r] << "," << count_alt[r] << ".";
				}

				topycc[rr] += sstr.str();
			}
		} // Get top (lowest level) group for the supposed derived allele
		top_groups.push_back(top_pair);

		// Work on coarser grouping
		for (int32_t l = 1; l <= mlevel; ++l) {
			std::map<std::string, int32_t> cref;
			std::map<std::string, int32_t> calt;
			for (int32_t i = 0; i < nbranch; ++i) {
				if (count_ref[i] > 0) {
					if (cref.find(branch_level[ branch[i] ][l]) == cref.end()) {
						cref.insert( { branch_level[ branch[i] ][l], count_ref[i] } );
					} else {
						cref[ branch_level[ branch[i] ][l] ] += count_ref[i];
					}
				}
				if (count_alt[i] > 0) {
					if (calt.find(branch_level[ branch[i] ][l]) == calt.end()) {
						calt.insert( { branch_level[ branch[i] ][l], count_alt[i] } );
					} else {
						calt[ branch_level[ branch[i] ][l] ] += count_alt[i];
					}
				}
			}
			ngr[l] = cref.size();
			nga[l] = calt.size();
			top_pair.clear();
			std::vector<std::pair<std::string, int32_t> > tpair;
			std::vector<std::pair<std::string, int32_t> > cvec;
			if (derived == 0) {
				for (auto & v : cref) {
					cvec.push_back({v.first, v.second});
				}
				// std::vector<std::pair<std::string, int32_t> > cvec = sort_map_val<std::string, int32_t>(cref);
			} else {
				for (auto & v : calt) {
					cvec.push_back({v.first, v.second});
				}
				// std::vector<std::pair<std::string, int32_t> > cvec = sort_map_val<std::string, int32_t>(calt);
			} // Get top group for the supposed derived allele at each level

			while(cvec.size() < 2) {
				cvec.push_back({"-", 0});
			}
			std::vector<int32_t> top_index(cvec.size());
			std::iota(top_index.begin(), top_index.end(), 0);
			sort(top_index.begin(), top_index.end(),
				[&cvec](int32_t left, int32_t right) -> bool {
						return cvec[left].second > cvec[right].second; } );
			for (int32_t rr = 0; rr < 2; ++rr) {
				int32_t r = top_index[rr];
				tpair.push_back({cvec[r]});
				std::stringstream sstr;
				sstr << cvec[r].first << "," << cvec[r].second << ".";
				topycc[rr] += sstr.str();
			}
			top_groups.push_back(tpair);
		}
		bcf_update_info_int32(odw.hdr, nv, "n_group_REF", &(ngr[0]), mlevel+1);
		bcf_update_info_int32(odw.hdr, nv, "n_group_ALT", &(nga[0]), mlevel+1);
		bcf_update_info_string(odw.hdr, nv, "YCCminor_1st", topycc[0].c_str());
		bcf_update_info_string(odw.hdr, nv, "YCCminor_2nd", topycc[1].c_str());
		odw.write(nv);
		bcf_destroy(nv);

		nVariant++;
		// if (ngr[olevel] <= 1 && nga[olevel] <= 1) {
		//   continue;
		// }
		std::string ref(1, iv->d.allele[0][0]);
		std::string alt(1, iv->d.allele[1][0]);
		std::string qual = ".";
		if (iv->qual == iv->qual) {
			qual = std::to_string((int) iv->qual);
		}
		wf << iv->pos << '\t' << ref << '\t' << alt << '\t' << qual << '\t';
		if (iv->d.n_flt == 0) {
			wf << '.';
		}
		else {
			wf << std::string( bcf_hdr_int2id(odr.hdr, BCF_DT_ID, iv->d.flt[0]) );
			if (iv->d.n_flt > 1) {
				for(int32_t i=1; i < iv->d.n_flt; ++i) {
					std::string filter_str( bcf_hdr_int2id(odr.hdr, BCF_DT_ID, iv->d.flt[i]) );
					wf << ',' << std::string( bcf_hdr_int2id(odr.hdr, BCF_DT_ID, iv->d.flt[i]) );
				}
			}
		}
		for (int32_t i = 0; i < nbranch; ++i) {
			wf << '\t' << count_ref[i] << '/' << count_alt[i];
		}
		wf << '\n';
	}

	wf.close();
	odw.close();
	return 0;
}
