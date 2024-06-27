#include "utils.h"
#include "cramore.h"
#include "hts_utils.h"
#include <chrono>
#include <random>
#include <numeric>
#include "haploid_hash.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"


int32_t cmdVcfCollapseVar(int32_t argc, char** argv) {

    std::string inVcf, samples, reg, out;
    int32_t verbose = 50000;
    bool haploid = false, diploid_homo = false, an_doubled = false;
    int32_t debug = 0;
    bool fold = false;
    bool preserve_order = false;
    bool sort_ac_decreasing = false;
    bool output_indiv_variant_info = false;
    int32_t min_ac = 2, max_ac = INT_MAX; // apply to minor ac
    int32_t sparse_thres = 10;
    double max_missing_rate = 0.1;
    bcf_vfilter_arg vfilt;
    // unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // std::mt19937 rng(seed);

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
        LONG_PARAM("snp-only",&vfilt.snpOnly,"Exclude indels")
        LONG_PARAM("diploid-homo",&diploid_homo,"Input is haploid though vcf is coded as diploid")
        LONG_PARAM("haploid",&haploid,"Input is haploid and the vcf is properly coded")
        LONG_PARAM("an-doubled",&an_doubled,"Input is haploid but AC/AN in the info field are doubled as the vcf is coded as diploid")
        LONG_PARAM("fold",&fold,"Flip the alleles so that variants with exactly complementary carriers are considered equal")
        LONG_INT_PARAM("min-ac",&min_ac,"Minimum minor allele count to consider")
        LONG_INT_PARAM("max-ac",&max_ac,"Maximum minor allele count to consider")
        LONG_INT_PARAM("sparse-thres",&sparse_thres,"Threshold for sparse representation")
        LONG_DOUBLE_PARAM("max-missing-rate",&max_missing_rate,"Maximum missing rate to keep a variant")

        LONG_PARAM_GROUP("Output Options", NULL)
        LONG_STRING_PARAM("out", &out, "Output file prefix")
        LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")
        LONG_PARAM("preserve-order",&preserve_order,"Preserve the order of the unique variants as the first appearance in the input VCF")
        LONG_PARAM("sort-ac-decreasing",&sort_ac_decreasing,"Sort the unique variants by decreasing ALT allele count. (Default is increasing)")
        LONG_PARAM("output-indiv-variant-info",&output_indiv_variant_info,"Output individual variant information")
        LONG_INT_PARAM("debug",&debug,"")
    END_LONG_PARAMS();
    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    // sanity check of input arguments
    if ( inVcf.empty() || out.empty()) {
    error("[E:%s:%d %s] --in-vcf, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
    }

    // Input VCF
    std::vector<GenomeInterval> intervals;
    if ( !reg.empty() ) {
        parse_intervals(intervals, "", reg);
    }
    BCFOrderedReader odr(inVcf, intervals);
    vfilt.init(odr.hdr);
    bcf1_t* iv = bcf_init();
    if (!samples.empty()) {
        bcf_hdr_set_samples(odr.hdr, samples.c_str(), 1);
    }
    int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
    uint32_t nalleles = 2 * nsamples;
    if (haploid || diploid_homo) {
        nalleles = nsamples;
    }
    int32_t pass_flt_id = bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "PASS");
    notice("Detected %d samples (%d haplotypes)", nsamples, nalleles);
    if (sparse_thres > nalleles / 64) {
        notice("Sparse threshold is too high. Reset to %d", nalleles / 64);
        sparse_thres = nalleles / 64;
    }

    // Store unique variants
    LshBinaryMix hmtx(nalleles, 32, sparse_thres, !fold);
    // Additional info
    std::vector<int32_t> obs_an;
    std::vector<int32_t> alt_ac;
    std::vector<int32_t> minor_ac; // incase we want to sort by minor allele
    std::vector<int32_t> npass;
    std::vector<bool>    flipped;
    std::vector<std::vector<int32_t> > positions;
    std::vector<std::vector<uint8_t> > filters; // Assuming there are at most 8 FILTER keys
    std::map<int32_t, int32_t> flt_idx;
    std::vector<std::string> flt_label;
    int32_t nflt = 0;
    for (uint32_t i = 0; i < odr.hdr->n[BCF_DT_ID]; ++i) {
        if (odr.hdr->id[BCF_DT_ID][i].key == NULL) continue;
        if (strcmp(odr.hdr->id[BCF_DT_ID][i].key, "PASS") == 0) {
            continue;
        }
        bcf_hrec_t* hrec = odr.hdr->id[BCF_DT_ID][i].val->hrec[BCF_HL_FLT];
        if (hrec == NULL) {continue;}
        flt_label.emplace_back(odr.hdr->id[BCF_DT_ID][i].key);
        for (uint32_t j = 0; j < hrec->nkeys; ++j) {
            if (strcmp(hrec->keys[j], "IDX") == 0) {
                flt_idx[i] = nflt++;
            }
        }
        if (nflt >= 8) {
            break;
        }
    }
    if (nflt > 8) {
        notice("Currently only support up to 8 FILTER keys for variant level output, %d keys are found.", nflt - 8);
    }

    int32_t  n_gt = 0, n_ac = 0, n_an = 0;
    int32_t* p_gt = NULL;
    int32_t* info_ac = NULL;
    int32_t* info_an = NULL;
    int32_t  k = 0, nVariant = 0, nUniq = 0, nPASS = 0, nFAIL = 0, nMiss = 0, nMissSkip = 0;
    int32_t nSparse = 0, nDense = 0;

    for(k=0; odr.read(iv); ++k) {  // read marker

        if ( k % verbose == 0 )
            notice("Processing %d records at %s:%d; %d contains missing gt, an additional %d contains more than %.1f %% missing gt and are skipped. Kept %d variants, recorded %d unique configurations. Query (hash) %d times, need to do full comparison %d times, %d times found exact match.", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nMiss, nMissSkip, max_missing_rate*100, nVariant, nUniq, hmtx.nQuery, hmtx.nCompare, hmtx.nCollision);
        bcf_unpack(iv, BCF_UN_SHR);
        if (!vfilt.pass_vfilter(iv)) {continue;}
        bool pass = false;
        for (size_t i = 0; i < iv->d.n_flt; ++i) {
            if (iv->d.flt[i] == pass_flt_id) {
                pass = true;
                break;
            }
        }
        if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
        if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
        if (an_doubled) {
            info_ac[0] /= 2;
            info_an[0] /= 2;
        }
        bool flip = info_ac[0] > info_an[0] / 2;
        int32_t ac = flip ? info_an[0] - info_ac[0] : info_ac[0]; // minor
        if (ac < min_ac || ac > max_ac) {continue;}
        bool sparse = ac <= sparse_thres;
        bcf_unpack(iv, BCF_UN_ALL);
        if (bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) <= 0 ) {
            error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
        }
        int32_t nmiss = 0, gt_ac = 0, gt_an = 0;
        int32_t ret = 0;
        int32_t j = diploid_homo ? 2 : 1;
        if (sparse) {
            std::set<uint32_t> carriers;
            for(int32_t i=0; i < nalleles; ++i) {
                if ( bcf_gt_is_missing(p_gt[i*j]) ) {
                    nmiss++;
                    continue;
                }
                gt_an++;
                gt_ac += bcf_gt_allele(p_gt[i*j]) > 0;
                if ((bcf_gt_allele(p_gt[i*j]) > 0) != flip) {
                    carriers.insert(i);
                }
            }
            if (nmiss > max_missing_rate * nalleles) {
                nMissSkip++;
                continue;
            }
            ac = carriers.size();
            if (ac < min_ac || ac > max_ac) {
                continue;
            }
            ret = hmtx.add_set(carriers, flip);
        } else {
            uint8_t gt[nalleles];
            for(int32_t i=0; i < nalleles; ++i) {
                if ( bcf_gt_is_missing(p_gt[i*j]) ) {
                    nmiss++;
                    gt[i] = flip ? 1 : 0;
                    continue;
                }
                gt_an++;
                if (bcf_gt_allele(p_gt[i*j]) > 0) {
                    gt[i] = 1;
                    gt_ac++;
                } else {
                    gt[i] = 0;
                }
            }
            if (nmiss > max_missing_rate * nalleles) {
                nMissSkip++;
                continue;
            }
            ac = gt_ac > gt_an / 2 ? gt_an - gt_ac : gt_ac;
            if (ac < min_ac || ac > max_ac) {
                continue;
            }
            ret = hmtx.add_int<uint8_t>(gt, false, flip);
        }
        uint8_t fltflag = 0;
        for(int32_t i=0; i < iv->d.n_flt; ++i) {

            if (flt_idx.find(iv->d.flt[i]) != flt_idx.end()) {
                fltflag |= 1 << flt_idx[iv->d.flt[i]];
            }
        }
        if (ret < 0) {
            ++nUniq;
            minor_ac.push_back(ac);
            alt_ac.push_back(gt_ac);
            obs_an.push_back(gt_an);
            npass.push_back(0);
            flipped.push_back(flip && (!fold));
            if (output_indiv_variant_info) {
                positions.push_back({(int32_t) iv->pos + 1});
                filters.push_back({fltflag});
            }

            ret = hmtx.nValues - 1;
            if (sparse) {
                nSparse++;
            } else {
                nDense++;
            }
        } else {
            if (output_indiv_variant_info) {
                positions[ret].push_back(iv->pos + 1);
                filters[ret].push_back(fltflag);
            }
        }
        if (nmiss > 0) {
            nMiss++;
        }
        if (!pass) {
            nFAIL++;
        } else {
            nPASS++;
            npass[ret]++;
        }
        ++nVariant;
if (debug % 3 == 1) {
std::cout << iv->pos+1 << " " << nUniq << "," << hmtx.nValues << '\t' << ret << '\t' << nmiss << "(" << nMissSkip << ")" << std::endl << std::flush;
}
        if (debug && nVariant > debug) {
            break;
        }
    }
    if (p_gt) free(p_gt);
    bcf_destroy(iv);
    odr.close();
    notice("Finished. Read %d records, kept %d variants, recorded %d unique configurations (%d sparse, %d dense). %d contains missing gt, an additional %d contains more than %.1f %% missing gt and are skipped. Query (hash) %d times, need to do full comparison %d times, %d times found exact match.", k, nVariant, nUniq, nSparse, nDense, nMiss, nMissSkip, max_missing_rate*100, hmtx.nQuery, hmtx.nCompare, hmtx.nCollision);

if (debug % 2 == 1) {
    return 0;
}
    std::vector<int32_t> indices(hmtx.nValues);
    std::iota(indices.begin(), indices.end(), 0);
    if (!preserve_order) {
        if (fold) {
            if (sort_ac_decreasing) {
                std::sort(indices.begin(), indices.end(), [&minor_ac](int i, int j) {return minor_ac[i] > minor_ac[j];});
            } else {
                std::sort(indices.begin(), indices.end(), [&minor_ac](int i, int j) {return minor_ac[i] < minor_ac[j];});
            }
        } else {
            if (sort_ac_decreasing) {
                std::sort(indices.begin(), indices.end(), [&alt_ac](int i, int j) {return alt_ac[i] > alt_ac[j];});
            } else {
                std::sort(indices.begin(), indices.end(), [&alt_ac](int i, int j) {return alt_ac[i] < alt_ac[j];});
            }
        }
    }

    // Output summary stats
    std::string outf = out + ".stats.tsv";
    FILE*  wf = fopen(outf.c_str(), "w");
    fprintf(wf, "AC\tAN\tnTotal\tnPASS\n");
    for (int32_t i = 0; i < hmtx.nValues; ++i) {
        int32_t idx = indices[i];
        fprintf(wf, "%d\t%d\t%d\t%d\n", alt_ac[idx], obs_an[idx], hmtx.nVariants[idx], npass[idx]);
    }
    fclose(wf);
    notice("Finished writing summary statistics.");

    if (output_indiv_variant_info) {
        outf = out + ".var_info.tsv.gz";
        htsFile* hf = hts_open(outf.c_str(), "wz");
        hprintf(hf, "POS\tConfigIndex\tFILTER_");
        for (auto it = flt_label.rbegin(); it != flt_label.rend(); ++it) {
            hprintf(hf, "%s|", it->c_str());
        }
        hprintf(hf, "\n");
        for (int32_t i = 0; i < hmtx.nValues; ++i) {
            int32_t idx = indices[i];
            for (size_t j = 0; j < positions[idx].size(); ++j) {
                hprintf(hf, "%d\t%d\t%d\n", positions[idx][j], i, filters[idx][j]);
            }
        }
        hts_close(hf);
    }

    // Output collapsed variant configurations
    outf = out + ".bin";
    wf = fopen(outf.c_str(), "wb");
    // meta data
    fwrite(&hmtx.nValues, sizeof(int32_t), 1, wf); // number of unique variants
    fwrite(&nalleles, sizeof(int32_t), 1, wf);     // number of haploid alleles
    fwrite(&sparse_thres, sizeof(int32_t), 1, wf); // threshold for sparse representation
    char delim[2] = {'>', '<'};                         // between variants
    char sparse_indicator[2] = {'0', '1'};         // 1 for sparse
    for (int32_t i = 0; i < hmtx.nValues; ++i) {
        int32_t idx = indices[i];
        if (flipped[idx]) {
            fwrite(delim+1, 1, 1, wf);                           // delimiter
        } else {
            fwrite(delim, 1, 1, wf);
        }
        fwrite(&minor_ac[idx], sizeof(int32_t), 1, wf);  // allele count
        fwrite(&hmtx.nVariants[idx], sizeof(int32_t), 1, wf); // total number
        fwrite(&npass[idx], sizeof(int32_t), 1, wf);          // number of PASS
        if (!hmtx.mixStorage[idx].isSparse) {
            fwrite(sparse_indicator, 1, 1, wf);               // 0 for dense
            fwrite(hmtx.mixStorage[idx].bitArray, 1, hmtx.nBytes, wf);
        } else {
            fwrite(sparse_indicator + 1, 1, 1, wf);           // 1 for sparse
            for (auto& v : hmtx.mixStorage[idx].indices) {
                fwrite(&v, sizeof(int32_t), 1, wf);
            }
        }
    }
    fclose(wf);
    notice("Finished dumpping full information.");

    return 0;
}
