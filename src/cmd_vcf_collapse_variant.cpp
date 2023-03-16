#include "utils.h"
#include "cramore.h"
#include "hts_utils.h"

#include <fstream>

#include "compact_matrix.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"

int32_t cmdVcfCollapseVar(int32_t argc, char** argv) {

  std::string inVcf, reg, out, samples, info;
  int32_t verbose = 10000;
  double minDprime = -1, maxMiss = 0.1;
  bool haploid = false;
  bool count_only = false;
  bool snp_only = false;
  bool anno_group = false;
  bool debug = false;
  double minAF = -1, maxAF = -1;
  int32_t minAC = 1, maxAC = -1;

  bcf_vfilter_arg vfilt;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
	LONG_PARAM_GROUP("Input Sites", NULL)
	LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
	LONG_STRING_PARAM("samples",&samples,"List of sample ID")
    LONG_STRING_PARAM("sample-info",&info,"Sample population or other group info")

	LONG_PARAM_GROUP("Variant Filtering Options", NULL)
	LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
	LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
	LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")
    LONG_PARAM("snp-only",&snp_only,"Exclude indels")
    LONG_DOUBLE_PARAM("min-af",&minAF,"Minimum minor allele frequency")
    LONG_DOUBLE_PARAM("max-af",&maxAF,"Maximum minor allele frequency")
    LONG_INT_PARAM("min-ac",&minAC,"Minimum minor allele count")
    LONG_INT_PARAM("max-ac",&maxAC,"Maximum minor allele count")
    LONG_PARAM("haploid",&haploid,"Input is haploid though vcf is read as diploid (?)")
    LONG_DOUBLE_PARAM("max-missing",&maxMiss,"Maximum missingness")

	LONG_PARAM_GROUP("Output Options", NULL)
	LONG_STRING_PARAM("out", &out, "Output file prefix")
	LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")
    LONG_PARAM("debug",&debug,"")



  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

    // sanity check of input arguments
    if ( inVcf.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
    }

    struct VarConfig {
        int32_t ac, n;
        std::map<std::string, int32_t> group_ct;
        std::map<int32_t, int32_t> filter_ct;
        std::string first_occur;
        std::string id;
    };

    // Output collapsed variant configurations
    std::string outf = out + ".var.config.tsv.gz";
    htsFile* wf = hts_open(outf.c_str(), "wg");
    std::map<int32_t, std::map<std::string, VarConfig> > clpsVars;

    std::vector<GenomeInterval> intervals;
    if ( !reg.empty() ) {
        parse_intervals(intervals, "", reg);
    }
    BCFOrderedReader odr(inVcf, intervals);
    bcf1_t* iv = bcf_init();
    int32_t pass_id = (int32_t) bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "PASS");

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
        std::string line;
        std::vector<std::string> words;
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
    if (minAF >= 0) {
        minAC = (int32_t) minAF * nalleles;
    }
    if (maxAF >= 0 && maxAF <= 0.5) {
        maxAC = (int32_t) maxAF * nalleles;
    }
    if (minAC < 0) { minAC = 1; }
    if (maxAC < 0) { maxAC = (int32_t) (nalleles / 2); }

    int32_t  n_gt = 0, n_ac = 0, n_an = 0;
    int32_t* p_gt = NULL;
    int32_t* info_ac = NULL;
    int32_t* info_an = NULL;
    int32_t  nVariant = 0, nUniq = 0;

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

        if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
        }
        if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
        if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
        int32_t ac = info_ac[0];
        if (info_an[0] > nalleles && haploid) {ac = ac / 2;}
        bool flip = ac > nalleles / 2;
        if (flip) {ac = nalleles - ac;}
        if (ac < minAC || ac > maxAC) {continue;}

        ac = 0;
        std::string config_id;
        std::stringstream ss;
        std::map<std::string, int32_t> group_count;
        for (const auto & v : group_list) {
            group_count[v] = 0;
        }
        int32_t ct = 0, n_miss = 0;
        int32_t current_gt = 0;
        int32_t current_start = 0;
        if (haploid) {
            if (!bcf_gt_is_missing(p_gt[0])) {
                current_gt = (bcf_gt_allele(p_gt[0]) > 0) ? 1 : 0;
                n_miss++;
            }
            for(int32_t i=1; i < nsamples; ++i) {
                uint8_t geno = 0;
                if ( bcf_gt_is_missing(p_gt[i*2]) ) {
                    geno = 0;
                    n_miss++;
                } else {
                    geno = (bcf_gt_allele(p_gt[i*2]) > 0) ? 1 : 0;
                }
                if (flip) {
                    geno = 1 - geno;
                    ac  += geno;
                } else {
                    ac += geno;
                }

                if (geno == 1) {
                    if (current_gt == 0) {
                        current_start = i;
                    }
                    ct += 1;
                }
                if (geno == 0 && current_gt == 1) {
                    ss << current_start << "+" << ct << "R;";
                    ct = 0;
                }
                current_gt = geno;
                group_count[sample_info[i]] += geno;
            }
            if (current_gt == 1)
                ss << current_start << "+" << ct << "R;";
        } else {
            if (!bcf_gt_is_missing(p_gt[0]) && !bcf_gt_is_missing(p_gt[1])) {
                current_gt = ((bcf_gt_allele(p_gt[0]) > 0) ? 1 : 0) + ((bcf_gt_allele(p_gt[1]) > 0) ? 1 : 0);
                n_miss++;
            }
            for(int32_t i=1; i < nsamples; ++i) {
                uint8_t geno = 0;
                if ( bcf_gt_is_missing(p_gt[i*2]) | bcf_gt_is_missing(p_gt[i*2+1]) ) {
                    geno = 0;
                    n_miss++;
                } else {
                    geno = ((bcf_gt_allele(p_gt[i*2]) > 0) ? 1 : 0) + ((bcf_gt_allele(p_gt[i*2+1]) > 0) ? 1 : 0);
                }
                if (flip) {
                    geno = 2 - geno;
                    ac  += geno;
                } else {
                    ac += geno;
                }

                if (geno != current_gt) {
                    if (current_gt == 1) {
                        ss << current_start << "+" << ct << "H;";
                    } else if (current_gt == 2) {
                        ss << current_start << "+" << ct << "A;";
                    }
                    ct = 0;
                    if (geno > 0) {
                        current_start = i;
                        ct = 1;
                    }
                } else if (geno > 0) {
                    ct += 1;
                }
                current_gt = geno;
                group_count[sample_info[i]] += geno;
            }
            if (current_gt > 0) {
                ss << current_start << "+" << ct;
                if (current_gt == 1) {
                    ss << "H;";
                } else {
                    ss << "A;";
                }
            }
        }
        if (n_miss > nsamples * maxMiss) {continue;}
        if (ac < minAC || ac > maxAC) {continue;}
        std::string ckey = std::to_string(ac) + "_" + ss.str();
        auto it = clpsVars.find(ac);
        if (it == clpsVars.end()) {
            clpsVars[ac] = std::map<std::string, VarConfig> {};
        }
        it = clpsVars.find(ac);
        auto ptr = (it->second).find(ckey);
        if (ptr == (it->second).end()) {
            VarConfig new_var;
            new_var.ac = ac;
            new_var.n  = 1;
            if (anno_group)
                new_var.group_ct = group_count;
            for(int32_t i=0; i < iv->d.n_flt; ++i)
                new_var.filter_ct[iv->d.flt[i]] = 1;

            ss.str(std::string());
            ss << ac << "_" << nUniq;
            new_var.id = ss.str();
            config_id = new_var.id;

            ss.str(std::string());
            ss << iv->pos+1 << "_" <<  bcf_get_allele(iv)[0] << "_" << bcf_get_allele(iv)[1];
            new_var.first_occur = ss.str();

            (it->second)[ckey] = new_var;
            nUniq += 1;
        } else {
            config_id = (ptr->second).id;
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
