#include "cramore.h"
#include "utils.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "ibs0.h"
#include "rare_variant_config.h"
#include <iomanip>
#include <numeric>
#include <sstream>

void UpdateInfo_IBS0(bcf_hdr_t *hdr, RareVariant *rare) {
  bcf1_t *nv = rare->iv;
  int32_t bc[2];
  int32_t ac = rare->ac;
  bc[0] = rare->left_size;
  bc[1] = ac - bc[0];
  bcf_update_info_int32(hdr, nv, "AvgDist_bp", &(rare->AvgDist), 1);
  bcf_update_info_float(hdr, nv, "AvgDist_cM", &(rare->AvgDist_cm), 1);
  bcf_update_info_float(hdr, nv, "MedDist_cM", &(rare->MedDist_cm), 1);
  bcf_update_info_int32(hdr, nv, "BiCluster", &bc, 2);
  bcf_update_info_int32(hdr, nv, "CarrierID", &(rare->id_list[0]), ac);
  bcf_update_info_float(hdr, nv, "Matrix_cM", &(rare->sorted_cm[0]), ac*ac);
  bcf_update_info_int32(hdr, nv, "Matrix_bp", &(rare->sorted_pt[0]), ac*ac);
}

void OutputRecord(BCFOrderedWriter& odw, RareVariant* rare, htsFile* wf, bp2cmMap& pgmap, double min_bp ) {
    rare->Organize(pgmap);
    UpdateInfo_IBS0(odw.hdr, rare);
    bcf_write(odw.file, odw.hdr, rare->iv);
    int32_t ac = rare->ac;
    // Output pairwise IBS0 intervals
    for (int32_t i = 0; i < ac-1; ++i) {
      for (int32_t j = i+1; j < ac; ++j) {
        int32_t r = rare->sorted_pt[i*ac+j] + rare->pos;
        int32_t l = rare->pos - rare->sorted_pt[j*ac+i];
        if (r - l + 1 >= min_bp) {
          hprintf(wf, "%d\t%d\t%d\t%d\t%.5f\n", rare->id_list[i], rare->id_list[j], l, r, pgmap.bpinterval2cm(l, r));
        }
      }
    }
    delete rare;
}

void copyInfoFields(bcf_hdr_t* src_hdr, bcf1_t* src_rec, bcf_hdr_t* dest_hdr, bcf1_t* dest_rec, const std::vector<std::string>& info_tags) {
    for (const auto& tag : info_tags) {
        int tag_id = bcf_hdr_id2int(src_hdr, BCF_DT_ID, tag.c_str());
        if (bcf_hdr_idinfo_exists(src_hdr, BCF_HL_INFO, tag_id)) {
            int32_t* val = NULL;
            int num_val = 0;
            int ret = bcf_get_info_int32(src_hdr, src_rec, tag.c_str(), &val, &num_val);
            if (ret > 0) {
                bcf_update_info_int32(dest_hdr, dest_rec, tag.c_str(), val, num_val);
            }
            if (val) free(val);
        }
    }
}

// goal -- get all pairwise no-ibs0 regions covering shared
//         rare variants from the given region
int32_t AnnotateIBS0AroundRare_Small(int32_t argc, char** argv) {
  std::string in_vcf, in_map, chrom, reg, oreg;
  std::string out_vcf, out_pair;
  int32_t min_hom_gts = 1;
  int32_t verbose = 1000;
  int32_t min_variant = 1;
  int32_t max_rare_ac = 10;
  int32_t min_rare_ac = 2;
  int32_t start, end;
  int32_t leftover = 0;
  int32_t cst = -1, ced = -1;
  int32_t ck_len = 500000;
  int32_t bp_limit = 5000000;
  double output_min_ibs0 = 5e5;
  bool rm_info = false;
  bool snp_only = false;
  std::vector<std::string> kept_info = {"AC","AN"};
  std::vector<std::string> info2remove;

  bcf_vfilter_arg vfilt;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&in_vcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
    LONG_STRING_PARAM("map",&in_map, "Map file for genetic distance")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")
    LONG_PARAM("snp-only",&snp_only, "Only consider SNPs")
    LONG_INT_PARAM("centromere-st",&cst, "Start position of centromere")
    LONG_INT_PARAM("centromere-ed",&ced, "End position of centromere")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("min-hom",&min_hom_gts, "Minimum number of homozygous genotypes to be counted for IBS0")
    LONG_INT_PARAM("min-variant",&min_variant, "Minimum number of variants to present to have output file")
    LONG_INT_PARAM("left-over",&leftover, "Stop when only x pairs are left")
    LONG_INT_PARAM("bp-limit",&bp_limit, "The upper bound of IBS0 search in each direction (bp). 0 for no limit")
    LONG_INT_PARAM("chunk-length",&ck_len, "The length of window to store common (0/1/2) variants (bp)")

    LONG_INT_PARAM("max-rare",&max_rare_ac, "Maximal minor allele count to be considered as anchor for IBD")
    LONG_INT_PARAM("min-rare",&min_rare_ac, "Minimal minor allele count to be considered as anchor for IBD")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_PARAM("rm-info",&rm_info, "Remove all INFO fields from the output except for those specified by --keep-info-tag")
    LONG_MULTI_STRING_PARAM("keep-info-tag",&kept_info, "Info tags to carry over to the output file")
    LONG_STRING_PARAM("out-vcf", &out_vcf, "Output VCF file name")
    LONG_STRING_PARAM("out-pair", &out_pair, "Output VCF file name")
    LONG_DOUBLE_PARAM("output-min-ibs0", &output_min_ibs0, "Minimum length of IBS0 to output separately (bp)")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( in_vcf.empty() || reg.empty() || in_map.empty() || out_vcf.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --region, --map, --out-vcf are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // Genetic map
  bp2cmMap pgmap(in_map, " ", "", cst, ced);
  notice("Read map. min %d; max %d; cst %d; ced %d.", pgmap.minpos, pgmap.maxpos, pgmap.centromere_st, pgmap.centromere_ed);
  if (bp_limit <= 0) {
    bp_limit = pgmap.maxpos - pgmap.minpos;
  }

  // Region to process
  oreg = reg;
  std::vector<std::string> v;
  split(v, ":-", reg);
  chrom = v[0];
  if (!str2int32(v[1], start) || !str2int32(v[2], end)) {
    error("Invalid region.");
  }

  // bcf reader
  // TODO: currently assume a small enough region.
  std::vector<GenomeInterval> intervals;
  parse_intervals(intervals, "", reg);
  BCFOrderedReader odr(in_vcf, intervals);
  bcf1_t* iv = bcf_init();

  // Output pairwise IBS0 intervals
  htsFile *wf = hts_open(out_pair.c_str(),"w");
  if (wf == NULL) {
    error("Failed to open file %s for writing", out_pair.c_str());
  }

  // bcf writer (site only)
  BCFOrderedWriter odw(out_vcf.c_str(),0);
  bcf_hdr_t* hnull = bcf_hdr_subset(odr.hdr, 0, 0, 0);
  bcf_hdr_remove(hnull, BCF_HL_FMT, NULL);
  if (rm_info) {
    for (int i = 0; i < odr.hdr->nhrec; ++i) {
        bcf_hrec_t* hrec = odr.hdr->hrec[i];
        // Check if the record is an INFO field
        if (hrec->type != BCF_HL_INFO) { continue; }
        for (int j = 0; j < hrec->nkeys; ++j) {
          if (std::string(hrec->keys[j]) != "ID") { continue; }
          if (std::find(kept_info.begin(), kept_info.end(), std::string(hrec->vals[j]) ) != kept_info.end()) {
            continue;
          }
          info2remove.push_back(hrec->vals[j]);
          // bcf_hdr_remove(hnull, BCF_HL_INFO, hrec->vals[j]);
          break;
        }
    }
  }
  odw.set_hdr(hnull);
  char buffer[65536];
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "BiCluster" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=BiCluster,Number=1,Type=String,Description=\"Sizes of the two subtree\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "CarrierID" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=CarrierID,Number=.,Type=Integer,Description=\"Carriers of the rare allele, order matches the IBD matrix\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "AvgDist_bp" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=AvgDist_bp,Number=1,Type=Integer,Description=\"Average no-IBS0 length between two clusters, in bp\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "AvgDist_cM" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=AvgDist_cM,Number=1,Type=Float,Description=\"Average no-IBS0 length between two clusters, in cM\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "MedDist_cM" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=MedDist_cM,Number=1,Type=Float,Description=\"Median no-IBS0 length between two clusters, in cM\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "Matrix_cM" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=Matrix_cM,Number=.,Type=Float,Description=\"No-IBS0 lengths (cM) of individual pairs arranged in a planer order; upper and lower triangular correspond to down- and upstream respectively\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "Matrix_bp" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=Matrix_bp,Number=.,Type=Integer,Description=\"Distance bewteen pairwise IBS0 positions to focal point arranged in a planer order; upper and lower triangular correspond to down- and upstream respectively\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  odw.write_hdr();

  // handle filter string
  std::string filter_str;
  int32_t filter_logic = 0;
  if ( vfilt.include_expr.empty() ) {
    if ( !vfilt.exclude_expr.empty() ) {
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

  int32_t nVariant = 0, nFinished = 0;
  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
  if ( nsamples == 0 )
    error("FATAL ERROR: The VCF does not have any samples with genotypes");
  notice("Input VCF file contains %d samples", nsamples);

  // Initialize IBS0lookup object surrounding the input region
  IBS0lookup ibs0finder(in_vcf, reg, pgmap, ck_len, ck_len, 1);
  if (ibs0finder.start_que.size() < 1) {
    notice("Not enough variants in the given region. Stopped without output.");
    return 0;
  }

  // Rare variants shared by pair of individuals
  //               where ibs0 range is not complete
  // individual ID pair -> list of rare variants they share
  std::map< std::pair<int32_t, int32_t>, std::vector<RareVariant*> > idpair_l;
  std::map< std::pair<int32_t, int32_t>, std::vector<RareVariant*> > idpair_r;
  // List of all focal rare variants
  // POS (1-based) -> ptr to a RareVariant object
  std::map<int32_t, RareVariant*> snplist;
  int32_t* p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t *info_ac = NULL;
  int32_t *info_an = NULL;
  int32_t n_ac = 0, n_an = 0, pos = 0;

  // First pass.
  // Find indiv. carrying rare alleles.
  for(int32_t k=0; odr.read(iv); ++k) {
    // periodic message to user
    if ( k % verbose == 0 )
      notice("Processing %d markers at %s:%d. Recorded %d rare variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nVariant);

    // unpack FILTER column
    bcf_unpack(iv, BCF_UN_SHR);
    if (snp_only && !bcf_is_snp(iv)) { continue; }
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
    if (info_ac[0] < 2) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
    int32_t ac = 0, an = 0;
    bool flip = false;
    ac = info_ac[0]; an = info_an[0];
    if (ac > an/2) {
      ac = an - ac;
      flip = true;
    }
    if (ac < min_rare_ac || ac > max_rare_ac) {continue;}
    if (iv->pos+1 == pos) {
      continue; // Jump over tri-allelic sites
    }

    // Get carrier IDs
    std::vector<int32_t> carry;
    for(int32_t i=0; i < nsamples; ++i) {
      int32_t g1 = p_gt[2*i];
      int32_t g2 = p_gt[2*i+1];
      int32_t geno;
      if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
        continue;
      } else {
        geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
        if ((!flip && geno >= 1) || (flip && geno <= 1) ) { // TODO: here we treat hom and het carriers equally
          carry.push_back(i);
        }
      }
    }
    int32_t n_carry = carry.size();
    if (n_carry < 2) {continue;}

    // Creat a new RareVariant record
    pos = iv->pos+1;
    bcf1_t* nv = bcf_dup(iv);
    bcf_unpack(nv, BCF_UN_SHR);
    bcf_subset(odw.hdr, nv, 0, 0);
    if (rm_info) {
      // Remove all INFO fields
      for (auto& tag : info2remove) {
        bcf_update_info(odw.hdr, nv, tag.c_str(), NULL, 0, BCF_HT_INT);
      }
    }
    RareVariant* rare = new RareVariant(nv, carry);
    snplist[pos] = rare;
    for (int32_t i = 0; i < n_carry-1; ++i) {
      for (int32_t j = i+1; j < n_carry; ++j) {
        std::pair<int32_t, int32_t> kpair(carry[i],carry[j]);
        idpair_r[kpair].push_back(rare);
        idpair_l[kpair].push_back(rare);
      }
    }

//     // Previous version, will find ibs0 repeatedly if involved in multiple rare vairants
//     // Find pairwise ibs0 w/in a short limit
//     int32_t fin = 0;
//     for (int32_t i = 0; i < n_carry-1; ++i) {
//       for (int32_t j = i+1; j < n_carry; ++j) {
//         int32_t r = ibs0finder.FindIBS0(carry[i],carry[j],pos,0);
//         int32_t l = ibs0finder.FindIBS0(carry[i],carry[j],pos,1);
//         if (l > 0 && r > 0) {
//           fin = rare->Add(carry[i], carry[j], l, r);
//         } else {
//           std::pair<int32_t, int32_t> kpair(carry[i],carry[j]);
//           if (l > 0) {
//             fin = rare->AddHalf(carry[i], carry[j], l, 1);
//             idpair_r[kpair].push_back(rare);
//           } else if (r > 0) {
//             fin = rare->AddHalf(carry[i], carry[j], r, 2);
//             idpair_l[kpair].push_back(rare);
//           } else {
//             idpair_r[kpair].push_back(rare);
//             idpair_l[kpair].push_back(rare);
//           }
//         }
//       }
//     }
//     if (fin) {
//       // Build tree
//       rare->Organize(pgmap);
// // std::cout << "Organized" << std::endl;
//       // Output this variant
//       UpdateInfo_IBS0(odw.hdr, rare);
// // std::cout << "Edit info" << std::endl;
//       bcf_write(odw.file, odw.hdr, nv);
// // std::cout << "Write record" << std::endl;
//       delete rare;
//       rare = NULL;
//       nFinished++;
//     } else { // Need to look further
//       snplist[pos] = rare;
//     }
    nVariant++;
  }

  // Find right ibs0 end points
  auto kv = idpair_r.begin();
  while (kv != idpair_r.end()) {
    auto itr = kv->second.begin();
    while(itr != kv->second.end()) {
      // Start from the left most focal rare variant
      int32_t r = ibs0finder.FindIBS0(kv->first.first,kv->first.second,(*itr)->pos,0);
      if (r < 0) {
        break;
      }
      while(itr != kv->second.end() && (*itr)->pos <= r) {
        // Add the right end point to all relevant variants
        int32_t fin = (*itr)->AddHalf(kv->first.first,kv->first.second,r,2);
        itr = kv->second.erase(itr);
      }
    }
    if (kv->second.size() == 0) {
      kv = idpair_r.erase(kv);
    } else {
      kv++;
    }
  }
  // Find left ibs0 end points
  kv = idpair_l.begin();
  while (kv != idpair_l.end()) {
    auto itr = kv->second.end();
    while(kv->second.size() > 0) {
      // Start from the right most focal rare variant
      itr--;
      int32_t l = ibs0finder.FindIBS0(kv->first.first,kv->first.second,(*itr)->pos,1);
      if (l < 0) {
        break;
      }
      int32_t fin = (*itr)->AddHalf(kv->first.first,kv->first.second,l,1);
      itr = kv->second.erase(itr);
      while(kv->second.size() > 0) {
        itr--;
        if ((*itr)->pos < l) {
          itr++;
          break;
        }
        fin = (*itr)->AddHalf(kv->first.first,kv->first.second,l,1);
        itr = kv->second.erase(itr);
      }
    }
    if (kv->second.size() == 0) {
      kv = idpair_l.erase(kv);
    } else {
      kv++;
    }
  }
  // Remove finished rare variants
  auto itr = snplist.begin();
  while (itr != snplist.end()) {
    RareVariant* rare = itr->second;
    if (rare->IfDone()) {
      OutputRecord(odw, rare, wf, pgmap, output_min_ibs0);
      itr = snplist.erase(itr);
      nFinished++;
    } else {
      itr++;
    }
  }

  notice("Finished first pass. Processed %d rare variants across %d samples; finished %d.", nVariant, nsamples, nFinished);
  notice("%d pairs miss left ibs0; %d pairs miss right ibs0.", idpair_l.size(), idpair_r.size());
  if (nVariant < 1) {
    notice("Not enough rare variants in the given region");
    odw.close();
    return 0;
  }
  int32_t bound1 = ibs0finder.start_que[0];
  int32_t bound2 = ibs0finder.posvec_que.back()->back();
  // Look backward
  notice("Start looking backward");
  int32_t wed = bound2;
  int32_t wst = bound1;
  int32_t pct = 0, fin = 0;
  std::string wreg;
  while (idpair_l.size() > 0 && !ibs0finder.reached_leftend && wst > start - bp_limit) {
    pct = 0;
    wed = wst - 1;
    wst = std::max(wed - ck_len, 0);
    wreg = chrom + ":" + std::to_string(wst) + "-" + std::to_string(wed);
    int32_t ret = ibs0finder.Update_Fixed(wreg);
    if (ret == 0) {
      continue;
    }
    if (ret < 0) {
      break;
    }
    auto itr = idpair_l.begin();
    while (itr != idpair_l.end()) {
      // Iterate over pairs of individuals missing left ibs0 pt
      int32_t l = ibs0finder.FindIBS0((*itr).first.first,(*itr).first.second,wed,1);
      if (l > 0) {
        for (auto ptr : (*itr).second) {
          // Add this pos to all relevant variants
          fin = ptr->AddHalf((*itr).first.first,(*itr).first.second,l,1);
          if (fin) {
            snplist.erase(ptr->iv->pos+1);
            OutputRecord(odw, ptr, wf, pgmap, output_min_ibs0);
            nFinished++;
          }
        }
        itr = idpair_l.erase(itr);
        pct++;
      } else {
        itr++;
      }
    }
    notice("Looking backward. Found %d more end points in %s; %d pairs still miss left ibs0; %d rare variants left.", pct, wreg.c_str(), idpair_l.size(), nVariant-nFinished);
    if ((int32_t) snplist.size() < leftover) {break;}
  }
  notice("Finish backward search, %d pairs did not find left IBS0 points", idpair_l.size());
  if (idpair_l.size() > 0) {
    // Add either the lower bound or the left end of the chromosome arm to the remaining pairs
    int32_t l = ibs0finder.leftend;
    if (!ibs0finder.reached_leftend) {
      if (ibs0finder.start_que.size() > 0) {
        l = ibs0finder.start_que[0];
      } else {
        l = wst;
      }
    }
    auto itr = idpair_l.begin();
    while (itr != idpair_l.end()) {
      for (auto ptr : (*itr).second) {
        fin = ptr->AddHalf((*itr).first.first,(*itr).first.second,l,1);
        if (fin) {
          snplist.erase(ptr->iv->pos+1);
          OutputRecord(odw, ptr, wf, pgmap, output_min_ibs0);
          nFinished++;
        }
      }
      itr = idpair_l.erase(itr);
    }
  }

  // Look forward
  notice("Start looking forward");
  wed = bound2;
  wst = bound1;
  while (idpair_r.size() > 0 && !ibs0finder.reached_rightend && wed < end + bp_limit) {
    pct = 0;
    wst = wed + 1;
    wed = wst + ck_len;
    wreg = chrom + ":" + std::to_string(wst) + "-" + std::to_string(wed);
    int32_t ret = ibs0finder.Update_Fixed(wreg);
    if (ret == 0) {
      continue;
    }
    if (ret < 0) {
      break;
    }
    auto itr = idpair_r.begin();
    while (itr != idpair_r.end()) {
      // Iterate over pairs of individuals missing left ibs0 pt
      int32_t r = ibs0finder.FindIBS0((*itr).first.first,(*itr).first.second,wst,0);
      if (r > 0) {
        for (auto ptr : (*itr).second) {
          // Add this pos to all relevant variants
          fin = ptr->AddHalf((*itr).first.first,(*itr).first.second,r,2);
          if (fin) {
            snplist.erase(ptr->iv->pos+1);
            OutputRecord(odw, ptr, wf, pgmap, output_min_ibs0);
            nFinished++;
          }
        }
        itr = idpair_r.erase(itr);
        pct++;
      } else {
        itr++;
      }
    }
    notice("Looking forward. Found %d more end points in %s; %d pairs still miss right ibs0; %d rare variants left.", pct, wreg.c_str(), idpair_r.size(), nVariant-nFinished);
    if ((int32_t) snplist.size() < leftover) {break;}
  }
  notice("Finish forward search, %d pairs did not find right IBS0 points", idpair_r.size());
  if (idpair_r.size() > 0) {
    int32_t r = ibs0finder.rightend;
    if (!ibs0finder.reached_rightend) {
      if (ibs0finder.posvec_que.size() > 0) {
        r = ibs0finder.posvec_que.back()->back();
      } else {
        r = wed;
      }
    }
    auto itr = idpair_r.begin();
    while (itr != idpair_r.end()) {
      for (auto ptr : (*itr).second) {
        fin = ptr->AddHalf((*itr).first.first,(*itr).first.second,r,2);
        if (fin) {
          snplist.erase(ptr->iv->pos+1);
          OutputRecord(odw, ptr, wf, pgmap, output_min_ibs0);
          nFinished++;
        }
      }
      itr = idpair_r.erase(itr);
    }
  }

  // Output unfinished variants (never happens)
  if (snplist.size() > 0) {
    notice("There are %d rare variants left.", snplist.size());
    for (auto ptr : snplist)
      bcf_write(odw.file, odw.hdr, ptr.second->iv);
  }

  odw.close();
  hts_close(wf);
  return 0;
}
