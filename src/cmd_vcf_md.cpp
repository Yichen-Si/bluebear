#include <stdio.h>
#include <string>
#include <sstream>
#include <set>
#include "cramore.h"
#include "utils.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"

struct popair {
		int32_t paternal_id;
		std::vector<int32_t> child_id;
		// (INDEL,SNP)x(FAIL,PASS)
		std::map<int32_t, std::vector<int32_t> > md_count_10;
		std::map<int32_t, std::vector<int32_t> > md_count_01;
};

int32_t MDHaploidY(int32_t argc, char** argv) {
	std::string inVcf, inPed, inSex;
	std::string out, outf;
	std::string reg;
	bool snp_only = false;
	int32_t debug = 0;
	int32_t verbose = 100000;
	int32_t n_family = 0, n_uniq_indiv = 0;

	bcf_vfilter_arg vfilt;

	paramList pl;
	BEGIN_LONG_PARAMS(longParameters)
		LONG_PARAM_GROUP("Input Sites", NULL)
		LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file of genotype")
		LONG_STRING_PARAM("in-ped",&inPed, "Input PED file of family info")
		LONG_STRING_PARAM("in-sex",&inSex, "Input file of sex info")
		LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

		LONG_PARAM_GROUP("Variant Filtering Options", NULL)
		LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
		LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
		LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")
		LONG_PARAM("snp-only",&snp_only,"Exclude indels")
		LONG_INT_PARAM("debug",&debug,"Debug")

		LONG_PARAM_GROUP("Output Options", NULL)
		LONG_STRING_PARAM("out", &out, "Output file prefix")
		LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

	END_LONG_PARAMS();

	pl.Add(new longParams("Available Options", longParameters));
	pl.Read(argc, argv);
	pl.Status();

	// sanity check of input arguments
	if ( inVcf.empty() || out.empty() || inPed.empty() ) {
		error("[E:%s:%d %s] --in-vcf, --in-ped and --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
	}

	// bcf reader
	std::vector<GenomeInterval> intervals;
	if ( !reg.empty() ) {
		parse_intervals(intervals, "", reg);
	}
	BCFOrderedReader odr(inVcf, intervals);
	bcf1_t* iv = bcf_init();
	int32_t pass_id = (int32_t) bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "PASS");

	int32_t n_samples = bcf_hdr_nsamples(odr.hdr);
	if ( n_samples == 0 )
		error("FATAL ERROR: The VCF does not have any samples with genotypes");
	char** id_ptr = bcf_hdr_get_samples(odr.hdr);
	std::map<std::string, int32_t> id_map;
	for (int32_t i = 0; i < n_samples; ++i) {
			id_map[std::string(id_ptr[i])] = i;
	}

	std::ifstream rf;
	std::string line;
	std::vector<std::string> v, u;
	// sex info
	std::set<std::string> male_id;
	rf.open(inSex, std::ifstream::in);
	while (std::getline(rf, line)) {
		split(v, "\t", line);
		if (v.size() < 2) {
			error("Invalid input sed information.");
		}
		auto p_it = id_map.find(v[0]);
		if (p_it == id_map.end()) {
			continue;
		}
		if (std::stoi(v[1]) == 1) {
			male_id.insert(v[0]);
		}
		v.clear();
	}
	rf.close();
	notice("Read %d males from --in-sex", male_id.size());

	// family info
	std::map<int32_t, popair> trios; // Paternal ID -> father-son pairs
	// std::vector<popair> trios;
	std::set<int32_t> all_po_id;
	rf.open(inPed, std::ifstream::in);
	while (std::getline(rf, line)) {
		split(v, "\t", line);
		if (v.size() < 4) {
			error("Invalid input family information.");
		}
		auto p_it = id_map.find(v[2]);
		if (p_it == id_map.end()) {
			continue;
		}
		popair rec = {p_it->second};
		split(u, ",", v[1]);
		for (const auto& w : u) { // Multiple IDs in the child column
			if (male_id.find(w) == male_id.end()) {
				continue; // Check if offspring is male
			}
			auto c_it = id_map.find(w);
			if (c_it == id_map.end()) {
				continue;
			}
			rec.child_id.push_back(c_it->second);
		}
		if (rec.child_id.size() == 0) {
			continue;
		}
		auto ptr_po = trios.find(p_it->second);
		if (ptr_po != trios.end()) { // This father is seem before
			for (const auto& w : rec.child_id) {
				all_po_id.insert(w);
				ptr_po->second.child_id.push_back(w);
			}
		} else {
			all_po_id.insert(rec.paternal_id);
			for (const auto& w : rec.child_id) {
				all_po_id.insert(w);
			}
			trios[p_it->second] = rec;
			n_family++;
		}

		if (debug and trios[p_it->second].child_id.size() > 1) {
			std::cout << n_family << '\t' << id_ptr[trios[p_it->second].paternal_id] << '\t' << trios[p_it->second].child_id.size() << '\n';
		}
		v.clear();
		u.clear();
	}
	rf.close();
	if ( n_family == 0 )
		error("FATAL ERROR: Input does not contain any trio");
	n_uniq_indiv = (int32_t) all_po_id.size();
	notice("Identified %d samples and %d families (including %d unique males as fathers or sons)", n_samples, n_family, n_uniq_indiv);

	if (debug == 3) {
		return 0;
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
	if ( bcf_hdr_id2int(odw.hdr, BCF_DT_ID, "nMD01") < 0 ) {
			sprintf(buffer,"##INFO=<ID=nMD01,Number=1,Type=Integer,Description=\"Number of father son pairs showing as (0, 1)\">\n");
			bcf_hdr_append(odw.hdr, buffer);
	}
	if ( bcf_hdr_id2int(odw.hdr, BCF_DT_ID, "nMD10") < 0 ) {
			sprintf(buffer,"##INFO=<ID=nMD10,Number=1,Type=Integer,Description=\"Number of father son pairs showing as (1, 0)\">\n");
			bcf_hdr_append(odw.hdr, buffer);
	}
	if ( bcf_hdr_id2int(odw.hdr, BCF_DT_ID, "nOffspring") < 0 ) {
			sprintf(buffer,"##INFO=<ID=nOffspring,Number=1,Type=Integer,Description=\"Number of offsprings carring the minor allele\">\n");
			bcf_hdr_append(odw.hdr, buffer);
	}
	if ( bcf_hdr_id2int(odw.hdr, BCF_DT_ID, "nPaternal") < 0 ) {
			sprintf(buffer,"##INFO=<ID=nPaternal,Number=1,Type=Integer,Description=\"Number of fathers carring the minor allele\">\n");
			bcf_hdr_append(odw.hdr, buffer);
	}
	if ( bcf_hdr_id2int(odw.hdr, BCF_DT_ID, "unrelMAC") < 0 ) {
			sprintf(buffer,"##INFO=<ID=unrelMAC,Number=1,Type=Integer,Description=\"Number of minor alleles excluding individuals from father son pairs\">\n");
			bcf_hdr_append(odw.hdr, buffer);
	}
	odw.write_hdr();


	// handle filter string
	vfilt.init(odr.hdr);
	vfilt.snpOnly = snp_only;

	int32_t *p_gt = NULL, *p_ac = NULL, *p_an = NULL;
	int32_t  n_gt = 0, n_ac = 0, n_an = 0, ac = 0, an = 0;
	int32_t  n_variant = 0;
	std::vector<int32_t> n_md_total(4, 0);

	for(int32_t k=0; odr.read(iv); ++k) {

		if ( k % verbose == 0 )
			notice("Processing vcf: %d markers at %s:%d. Found %d mendelian discordance cases among PASS variants (%d SNP, %d INDEL)", n_variant, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, n_md_total[1]+n_md_total[3], n_md_total[3], n_md_total[1]);
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

		bool is_pass = false;
		for(int32_t i=0; i < iv->d.n_flt; ++i) {
			if (iv->d.flt[i] == pass_id)
				is_pass = true;
				break;
		}
		int32_t variant_type = (int32_t) bcf_is_snp(iv) * 2 + (int32_t) is_pass;

		if (bcf_get_info_int32(odr.hdr, iv, "AC", &p_ac, &n_ac) < 0) {continue;}
		if (bcf_get_info_int32(odr.hdr, iv, "AN", &p_an, &n_an) < 0) {continue;}
		ac = p_ac[0];
		an = p_an[0];
		int32_t minor = 1;
		if (ac == 0 || ac == an) { continue; }
		if (ac > an / 2) {
			minor = 0;
		}
		if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
			error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
		}
		bool haploid = 0;
		if (n_gt == n_samples) {
			haploid = 1;
		}
		if (debug == 4) {
			std::cout << n_gt << '\t' << haploid << '\t' << n_samples << '\n';
			return 0;
		}
		int32_t ac_exclude_family = 0;
		int32_t ac_minor = 0;
		for (int32_t i = 0; i < n_samples; ++i) {
			int32_t gt = 0;
			if (haploid) {
				if ( bcf_gt_is_missing(p_gt[i]) ) {
					continue;
				} else {
					gt = (bcf_gt_allele(p_gt[i]) > 0) ? 1 : 0;
				}
			} else {
				if ( bcf_gt_is_missing(p_gt[2*i]) ) {
					continue;
				} else {
					gt = (bcf_gt_allele(p_gt[2*i]) > 0) ? 1 : 0;
				}
			}
			if (gt == minor) {
				ac_minor++;
				if (all_po_id.find(i) == all_po_id.end()) {
					ac_exclude_family++;
				}
			}
		}

		int32_t n_minor_offspring = 0, n_minor_paternal = 0;
		int32_t n_md_01 = 0, n_md_10 = 0;  // Number MD pairs
		for (auto& ptr : trios) {
			popair& rec = ptr.second;
			int32_t paternal  = -1;

			if (haploid) {
				if ( bcf_gt_is_missing(p_gt[rec.paternal_id]) ) {
					continue;
				} else {
					paternal = (bcf_gt_allele(p_gt[rec.paternal_id]) > 0) ? 1 : 0;
				}
			} else {
				if ( bcf_gt_is_missing(p_gt[2*rec.paternal_id]) ) {
					continue;
				} else {
					paternal = (bcf_gt_allele(p_gt[2*rec.paternal_id]) > 0) ? 1 : 0;
				}
			}

			int32_t md_10 = 0; // Count 10 (0 if father=0 or all pairs are 11)
			int32_t md_01 = 0; // Count 01 (0 if father=1 or all pairs are 00)
			n_minor_paternal += (int32_t) (paternal == minor);

			for (const auto& w : rec.child_id) {
				int32_t offspring = -1;
				if (haploid) {
					if ( bcf_gt_is_missing(p_gt[w]) ) {
						continue;
					} else {
						offspring = (bcf_gt_allele(p_gt[w]) > 0) ? 1 : 0;
					}
				} else {
					if ( bcf_gt_is_missing(p_gt[2*w]) ) {
						continue;
					} else {
						offspring = (bcf_gt_allele(p_gt[2*w]) > 0) ? 1 : 0;
					}
				}
				if (paternal == minor && offspring != minor) {
					md_10 = 1;
				}
				if (paternal != minor && offspring == minor) {
					md_01++;
					auto ptr = rec.md_count_01.find(ac_minor);
					if (ptr == rec.md_count_01.end()) {
						std::vector<int32_t> v(4, 0);
						v[variant_type] = 1;
						rec.md_count_01[ac_minor] = v;
					} else {
						ptr->second[variant_type]++;
					}
				}
				if (offspring == minor) {
					n_minor_offspring++;
				}

if ((debug+1) % 5 == 0 && md_10 + md_01 > 0) {
std::cout << paternal << '\t' << offspring << '\t' << ac_minor << '\t' << md_10 << ',' << md_01 << '\n';
}
if ((debug+1) % 6 == 0 && ac_minor == 1 && n_minor_paternal > 0 && paternal == minor && md_10 == 0) {
std::cout << iv->pos+1 << '\t' << n_minor_paternal << '\t' << paternal << '\t' << offspring << '\t' << ac_minor << '\t' << ac << '\t' << md_10 << ',' << md_01 << '\n';
}

			}
			if (md_10 > 0) {
				auto ptr = rec.md_count_10.find(ac_minor);
				if (ptr == rec.md_count_10.end()) {
					std::vector<int32_t> v(4, 0);
					v[variant_type] = 1;
					rec.md_count_10[ac_minor] = v;
				} else {
					ptr->second[variant_type]++;
				}
			}
			n_md_10 += md_10;
			n_md_01 += md_01;
		}
		n_md_total[variant_type] += n_md_10 + n_md_01;
		if ((debug == 3) && n_md_10 + n_md_01 > 0) {
			std::cout << ac_minor << '\t' << n_minor_offspring << '\t' << n_minor_paternal << '\t' << n_md_10 << ',' << n_md_01 << '\n';
		}

		// Output
		bcf1_t* nv = bcf_dup(iv);
		bcf_unpack(nv, BCF_UN_ALL);
		bcf_subset(odw.hdr, nv, 0, 0);
		bcf_update_info_int32(odw.hdr, nv, "nMD01", &n_md_01, 1);
		bcf_update_info_int32(odw.hdr, nv, "nMD10", &n_md_10, 1);
		bcf_update_info_int32(odw.hdr, nv, "nOffspring", &n_minor_offspring, 1);
		bcf_update_info_int32(odw.hdr, nv, "nPaternal", &n_minor_paternal, 1);
		bcf_update_info_int32(odw.hdr, nv, "unrelMAC", &ac_exclude_family, 1);
		odw.write(nv);
		bcf_destroy(nv);
		if (debug > 0 && n_md_total[3] > debug) {
			break;
		}
		n_variant++;
	}
	if (p_gt) free(p_gt);
	odw.close();

	// Output family summary
	outf = out + ".family_summary.tsv";
	FILE* wf = fopen(outf.c_str(), "w");
	fprintf(wf, "Paternal_id\tChild_id\tType\tMinor_AC\tnMD_PASS_SNP\tnMD_FAIL_SNP\tnMD_PASS_INDEL\tnMD_FAIL_INDEL\n");
	for (auto& po_ptr : trios) {
		popair& rec = po_ptr.second;
		std::stringstream ss;
		ss << id_ptr[rec.child_id[0]];
		if (rec.child_id.size() > 1) {
			for (uint32_t i = 1; i < rec.child_id.size(); ++i) {
				ss << ',' << id_ptr[rec.child_id[i]];
			}
		}
		for (const auto& ptr: rec.md_count_10) {
			fprintf(wf, "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\n",
			id_ptr[rec.paternal_id], ss.str().c_str(), "1/0", ptr.first,
			ptr.second[3], ptr.second[2], ptr.second[1], ptr.second[0]);
		}
		for (const auto& ptr: rec.md_count_01) {
			fprintf(wf, "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\n",
			id_ptr[rec.paternal_id], ss.str().c_str(), "0/1", ptr.first,
			ptr.second[3], ptr.second[2], ptr.second[1], ptr.second[0]);
		}
	}
	fclose (wf);
	odr.close();
	notice("Finish processing %d markers. Found %d mendelian discordance cases among PASS variants (%d SNP, %d INDEL)", n_variant, n_md_total[1]+n_md_total[3], n_md_total[3], n_md_total[1]);

	return 0;
}
