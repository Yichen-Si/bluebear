#include "cramore.h"
#include "utils.h"
#include "bcf_filter_arg.h"
#include "tsv_reader.h"
#include "ibs0.h"
#include "rare_variant_config.h"
#include <iomanip>
#include <numeric>
#include <sstream>
#include <iostream>

int32_t IBS0ApproxPMdistr(int32_t argc, char** argv) {

    std::string in_vcf, in_target, in_map, chrom, reg;
    std::string output;
    int32_t min_hom_gts = 1;
    int32_t verbose = 1000;
    int32_t start = -1, end = -1;
    int32_t wst, wed;
    int32_t cst = -1, ced = -1;
    int32_t ck_len = 500000;
    int32_t bp_limit = 2000000;
    int32_t window_size = 1000000;
    int32_t minpos = -1, maxpos = -1;
    int32_t minac = 20;

    paramList pl;

    BEGIN_LONG_PARAMS(longParameters)
        LONG_PARAM_GROUP("Input", NULL)
        LONG_STRING_PARAM("in-vcf",&in_vcf, "Input VCF/BCF file")
        LONG_STRING_PARAM("in-target",&in_target, "Input TSV file including paris (of sets) of sample IDs")
        LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
        LONG_STRING_PARAM("map",&in_map, "Map file for genetic distance")

        LONG_PARAM_GROUP("Variant Filtering Options", NULL)
        // LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
        // LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
        // LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")
        LONG_INT_PARAM("centromere-st",&cst, "Start position of centromere")
        LONG_INT_PARAM("centromere-ed",&ced, "End position of centromere")
        LONG_INT_PARAM("min-pos",&minpos, "Start position of the sequenced chromosome")
        LONG_INT_PARAM("max-pos",&maxpos, "End position of the sequenced chromosome")
        LONG_INT_PARAM("min-ac-hom",&minac, "Minimum AC to be considered for IBS0 calculation")

        LONG_PARAM_GROUP("Additional Options", NULL)
        LONG_INT_PARAM("min-hom",&min_hom_gts, "Minimum number of homozygous genotypes to be counted for IBS0")
        LONG_INT_PARAM("bp-limit",&bp_limit, "The upper bound of IBS0 search in each direction (bp). 0 for no limit")
        LONG_INT_PARAM("window-size",&window_size, "The size of the window to focus on at one time (bp)")
        LONG_INT_PARAM("chunk-length",&ck_len, "The length of window to store common (0/1/2) variants (bp)")

        LONG_PARAM_GROUP("Output Options", NULL)
        LONG_STRING_PARAM("out", &output, "Output file name")
        LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    // sanity check of input arguments
    if ( in_vcf.empty() || reg.empty() || in_map.empty() || output.empty() ) {
        error("[E:%s:%d %s] --in-vcf, --region, --map, --output are required parameters",__FILE__,__LINE__,__FUNCTION__);
    }

    // Genetic map
    bp2cmMap pgmap(in_map, " ", "", cst, ced);
    notice("Read map. min %d; max %d; cst %d; ced %d.", pgmap.minpos, pgmap.maxpos, pgmap.centromere_st, pgmap.centromere_ed);
    if (bp_limit <= 0) {
        bp_limit = pgmap.maxpos - pgmap.minpos;
    }
    if (minpos < 0) {
        minpos = pgmap.minpos;
    }
    if (maxpos < 0) {
        maxpos = pgmap.maxpos;
    }

    // Region to process
    std::vector<std::string> v;
    split(v, ":-", reg);
    chrom = v[0];
    if (v.size() > 1 && !str2int32(v[1], start)) {
        error("Invalid region.");
    }
    if (start < 0) {
        start = minpos;
    }
    if (v.size() > 2 && !str2int32(v[2], end)) {
        error("Invalid region.");
    }
    if (end < 0) {
        end = maxpos;
    }

    // Initialize target pairs
    tsv_reader tr(in_target.c_str());
    if ( tr.read_line() == 0 ) {
        error("[E:%s] Cannot read input query %s", __PRETTY_FUNCTION__,in_target.c_str());
    }

    // Read sample IDs from input VCF
    htsFile* fp = hts_open(in_vcf.c_str(), "r");
    if (fp == NULL) {
        error("Failed to open %s", in_vcf.c_str());
    }
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    int32_t nsamples = bcf_hdr_nsamples(hdr);
    std::map<std::string, int32_t> sample_idx;
    for (int32_t i = 0; i < nsamples; ++i) {
        sample_idx[bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i)] = i;
    }
    notice("Read %d sample IDs from %s", nsamples, in_vcf.c_str());

    // Output pairwise IBS0 intervals
    htsFile *wf = hts_open(output.c_str(),"wz");
    if (wf == NULL) {
        error("Failed to open %s for writing", output.c_str());
    }

    int32_t pos = 0, pos1 = 0, pos2 = 0;
    int32_t nrec = 0, nskip = 0;
    std::vector<std::string> id1, id2;
    std::vector<int32_t> idx1, idx2;
    std::vector<int32_t> bpvec;
    std::vector<double> cmvec;

    wst = start;
    // Initialize IBS0lookup object
    std::string reg1 = chrom + ":" + std::to_string(wst) + "-" + std::to_string(wst + window_size - 1);
    IBS0lookup ibs0finder(in_vcf, reg1, pgmap, bp_limit, ck_len, min_hom_gts, minpos, maxpos, minac);

    while (wst < end) {
        wed = wst + window_size - 1;
        if (wed > end) {
            wed = end;
        }
        reg1 = chrom + ":" + std::to_string(wst) + "-" + std::to_string(wed);
        bool ret = tr.jump_to(reg1.c_str());
        if (!ret) {
            notice("No query in region %s", reg1.c_str());
            wst = wed + 1;
            continue;
        }
        int32_t nret = ibs0finder.Update(reg1);
        if (nret < 0) {
            notice("Skip region %s", reg1.c_str());
            wst = wed + 1;
            continue;
        }
        notice("Processing region %s", reg1.c_str());

        while (tr.read_line() > 0) {
            idx1.clear();
            idx2.clear();
            bpvec.clear();
            cmvec.clear();
            if ((nrec + nskip) % verbose == 0) {
                notice("Processed %d queries. Skipped %d.", nrec, nskip);
            }
            pos1 = tr.int_field_at(1);
            pos2 = tr.int_field_at(2);
            pos = (pos1 + pos2) / 2;
            split(id1, ",", tr.str_field_at(3));
            split(id2, ",", tr.str_field_at(4));
            if (id1.size() < 1 || id2.size() < 1) {
                error("Incomplete IDs in the target file at line %s", tr.str.s);
            }
            bool flag = false;
            for (auto v : id1) {
                auto ptr = sample_idx.find(v);
                if (ptr == sample_idx.end()) {
                    continue;
                }
                idx1.push_back(ptr->second);
            }
            for (auto v : id2) {
                auto ptr = sample_idx.find(v);
                if (ptr == sample_idx.end()) {
                    continue;
                }
                idx2.push_back(ptr->second);
            }
            int32_t s1 = idx1.size(), s2 = idx2.size();
            if (s1 < 1 || s2 < 1) {
                nskip++;
                warning("At least one sample set is absent from the input vcf\ntarget line: %s", tr.str.s);
                continue;
            }
            for (auto i : idx1) {
                for (auto j : idx2) {
                    if (i == j) {
                        flag = true;
                        break;
                    }
                    int32_t r = ibs0finder.FindIBS0(i, j, pos, 0, 1);
                    int32_t l = ibs0finder.FindIBS0(i, j, pos, 1, 1);
                    if (r < 0 || l < 0) {
                        warning("Something is wrong - requested position is out of range (%d, %d, %d)", pos, l, r);
                        flag = true;
                        break;
                    }
                    bpvec.push_back(r - l);
                    cmvec.push_back(pgmap.bpinterval2cm(l, r));
                }
                if (flag) {
                    warning("Sample overlap (%s, %s)", tr.str_field_at(3), tr.str_field_at(4));
                    break;
                }
            }
            if (flag || bpvec.size() < 1) {
                nskip++;
                continue;
            }

            hprintf(wf, "%s\t%d\t%d\t%d\t", chrom.c_str(), pos, s1, s2);
            switch (bpvec.size()) {
                case 1:
                    hprintf(wf, "%d\t%.5f\t.\n", bpvec[0], cmvec[0]);
                    break;
                case 2:
                    hprintf(wf, "%d\t%.5f\t%d,%.5f\n", (bpvec[0]+bpvec[1])/2, (cmvec[0]+cmvec[1])/2, bpvec[0], cmvec[0]);
                    break;
                default:
                    std::sort(bpvec.begin(), bpvec.end());
                    std::sort(cmvec.begin(), cmvec.end());
                    int32_t midbp = bpvec[bpvec.size()/2];
                    double midcm = cmvec[cmvec.size()/2];
                    if (bpvec.size() % 2 == 0) {
                        midbp = (bpvec[bpvec.size()/2-1] + midbp) / 2;
                        midcm = (cmvec[cmvec.size()/2-1] + midcm) / 2;
                    }
                    hprintf(wf, "%d\t%.5f\t%d,%.5f,%d,%.5f,%d,%.5f\n",
                        midbp, midcm,
                        std::accumulate(bpvec.begin(), bpvec.end(), 0)  / bpvec.size(),
                        std::accumulate(cmvec.begin(), cmvec.end(), 0.) / cmvec.size(),
                        bpvec[0], cmvec[0], bpvec[bpvec.size()-1], cmvec[cmvec.size()-1]);
            }
            nrec++;
        }
        wst = wed + 1;
    }

    notice("Finished. Processed %d queries. Skipped %d.", nrec, nskip);
    hts_close(wf);
    return 0;
}
