#include "cramore.h"
#include "utils.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include <iomanip>
#include "ibs0.h"
#include "rare_variant_ibs0.h"

void UpdateInfo_IBS0(bcf_hdr_t *hdr, bcf1_t *nv, RareVariant *rare) {
  std::stringstream ss;
  std::string sstr;
  ss << rare->subset1.size() << "," << rare->subset2.size();
  bcf_update_info_string(hdr, nv, "BiCluster", ss.str().c_str());

  bcf_update_info_int32(hdr, nv, "LeftSet", &(rare->subset1[0]), rare->subset1.size());
  bcf_update_info_int32(hdr, nv, "RightSet", &(rare->subset2[0]), rare->subset2.size());
  bcf_update_info_int32(hdr, nv, "Carriers", &(rare->id_list[0]), rare->id_list.size());
  bcf_update_info_int32(hdr, nv, "AvgDist_bp", &(rare->AvgDist), 1);
  bcf_update_info_float(hdr, nv, "AvgDist_cM", &(rare->AvgDist_cm), 1);
  bcf_update_info_float(hdr, nv, "MedDist_cM", &(rare->MedDist_cm), 1);
  bcf_update_info_float(hdr, nv, "BtwDist", &(rare->btwDist), (rare->btwDist).size());

  ss.str(std::string());
  for (int32_t i = 0; i < rare->ac-1; ++i) {
    for (int32_t j = i+1; j < rare->ac; ++j) {
      ss << rare->ibs0mat[j][i] << ',';
    }
  }
  sstr = ss.str(); sstr.pop_back();
  bcf_update_info_string(hdr, nv, "LeftIBS0", sstr.c_str());

  ss.str(std::string());
  for (int32_t i = 0; i < rare->ac-1; ++i) {
    for (int32_t j = i+1; j < rare->ac; ++j) {
      ss << rare->ibs0mat[i][j] << ',';
    }
  }
  sstr = ss.str(); sstr.pop_back();
  bcf_update_info_string(hdr, nv, "RightIBS0", sstr.c_str());
}

// goal -- get all pairwise no-ibs0 regions covering shared
//         rare variants from the given region
int32_t AnnotateIBS0AroundRare_Small(int32_t argc, char** argv) {
  std::string inVcf, inMap, chrom, reg, oreg;
  std::string out;
  int32_t min_hom_gts = 1;
  int32_t verbose = 1000;
  int32_t min_variant = 1;
  int32_t max_rare_ac = 10;
  int32_t start, end;
  int32_t leftover = 0;
  int32_t cst = -1, ced = -1;
  int32_t ck_len = 500000;
  int32_t bp_limit = 1000000;
  // double  cm_limit =2.0;

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")
    LONG_STRING_PARAM("map",&inMap, "Map file for genetic distance")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")
    LONG_INT_PARAM("centromere-st",&cst, "Start position of centromere")
    LONG_INT_PARAM("centromere-ed",&ced, "End position of centromere")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("min-hom",&min_hom_gts, "Minimum number of homozygous genotypes to be counted for IBS0")
    LONG_INT_PARAM("min-variant",&min_variant, "Minimum number of variants to present to have output file")
    LONG_INT_PARAM("left-over",&leftover, "Stop when only x pairs are left")
    LONG_INT_PARAM("bp-limit",&bp_limit, "The length of window to store common (0/1/2) variants (bp)")

    // LONG_DOUBLE_PARAM("out-reach-cm",&cm_limit, "How far to look forward & backward (cm)")

    LONG_INT_PARAM("max-rare",&max_rare_ac, "Maximal minor allele count to be considered as anchor for IBD")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || reg.empty() || inMap.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --region, --inMap, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // Genetic map
  bp2cmMap pgmap(inMap, " ", "", cst, ced);
  notice("Read map. min %d; max %d; cst %d; ced %d.", pgmap.minpos, pgmap.maxpos, pgmap.centromere_st, pgmap.centromere_ed);

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
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();

  // bcf writer (site only)
  BCFOrderedWriter odw(out.c_str(),0);
  bcf_hdr_t* hnull = bcf_hdr_subset(odr.hdr, 0, 0, 0);
  bcf_hdr_remove(hnull, BCF_HL_FMT, NULL);
  odw.set_hdr(hnull);
  char buffer[65536];
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "BiCluster" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=BiCluster,Number=1,Type=String,Description=\"Sizes of the two subtree\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "LeftSet" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=LeftSet,Number=.,Type=Integer,Description=\"Individual ID in one subset\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "RightSet" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=RightSet,Number=.,Type=Integer,Description=\"Individual ID in another subset\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "Carriers" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=Carriers,Number=.,Type=Integer,Description=\"Individual ID carrying the rare allele\">\n");
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
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "LeftIBS0" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=LeftIBS0,Number=1,Type=String,Description=\"Pairwise IBS0 position to the left\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "RightIBS0" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=RightIBS0,Number=1,Type=String,Description=\"Pairwise IBS0 position to the right\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  if ( bcf_hdr_id2int(odr.hdr, BCF_DT_ID, "BtwDist" ) < 0 ) {
    sprintf(buffer,"##INFO=<ID=BtwDist,Number=.,Type=Float,Description=\"No-IBS0 lengths of individual pairs in the two final clusters, in cM\">\n");
    bcf_hdr_append(odw.hdr, buffer);
  }
  odw.write_hdr();

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

  int32_t nVariant = 0, nFinished = 0;
  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);

  notice("Started Reading site information from VCF file, identifying %d samples", nsamples);
  if ( nsamples == 0 )
    error("FATAL ERROR: The VCF does not have any samples with genotypes");

  IBS0lookup ibs0finder(inVcf, reg, pgmap, bp_limit/2, ck_len, 1);
  if (ibs0finder.start_que.size() < 1) {
    notice("Not enough variants in the given region. Stopped without output.");
    return 0;
  }

  std::map< std::pair<int32_t, int32_t>, std::vector<RareVariant*> > idpair_l;
  std::map< std::pair<int32_t, int32_t>, std::vector<RareVariant*> > idpair_r;
  std::map<int32_t, RareVariant*> snplist;
  int32_t* p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t *info_ac = NULL;
  int32_t *info_an = NULL;
  int32_t n_ac = 0, n_an = 0, pos = 0;

  // First pass.
  // Find indiv. carrying rare alleles.
  for(int32_t k=0; odr.read(iv); ++k) {
    if (iv->pos+1 == pos) {
      continue; // Jump over tri-allelic sites
    }
    // periodic message to user
    if ( k % verbose == 0 )
      notice("Processing %d markers at %s:%d. Recorded %d rare variants", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nVariant);

    // unpack FILTER column
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
    // extract genotype and apply genotype level filter
    if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (info_ac[0] < 2) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
    int32_t ac = 0, an = 0;
    ac = info_ac[0]; an = info_an[0];
    if (ac > an/2) {
      ac = an - ac;
    }
    if (ac < 2 || ac > max_rare_ac) {continue;}

    std::vector<int32_t> carry;
    for(int32_t i=0; i < nsamples; ++i) {
      int32_t g1 = p_gt[2*i];
      int32_t g2 = p_gt[2*i+1];
      int32_t geno;
      if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
        continue;
      } else {
        geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
        if (geno == 1) {
          carry.push_back(i);
        }
      }
    } // Found carriers (have to be het)
    // TODO: here we ignore rare varaiants with any hom carrier
    if ((int32_t) carry.size() != ac) {continue;}

    // Creat a new RareVariant record
    pos = iv->pos+1;
    bcf1_t* nv = bcf_dup(iv);
    bcf_unpack(nv, BCF_UN_ALL);
    bcf_subset(odw.hdr, nv, 0, 0);
    RareVariant* rare = new RareVariant(nv, ac);
    rare->Add_id(carry);
    // Find pairwise ibs0 w/in a short limit
    int32_t fin = 0;
    for (int32_t i = 0; i < ac-1; ++i) {
      for (int32_t j = i+1; j < ac; ++j) {
        int32_t r = ibs0finder.FindIBS0(carry[i],carry[j],pos,0);
        int32_t l = ibs0finder.FindIBS0(carry[i],carry[j],pos,1);
        if (l > 0 && r > 0) {
          fin = rare->Add(carry[i], carry[j], l, r);
          continue;
        } else {
          std::pair<int32_t, int32_t> kpair(carry[i],carry[j]);
          if (l > 0) {
            fin = rare->Add_half(carry[i], carry[j], l, 1);
            idpair_r[kpair].push_back(rare);
          } else if (r > 0) {
            fin = rare->Add_half(carry[i], carry[j], r, 2);
            idpair_l[kpair].push_back(rare);
          } else {
            idpair_r[kpair].push_back(rare);
            idpair_l[kpair].push_back(rare);
          }
        }
      }
    }
    if (fin) {
      // Build tree
      rare->Organize(pgmap);
      // Output this variant
      UpdateInfo_IBS0(odw.hdr, nv, rare);
      odw.write(nv);
      delete rare;
      rare = NULL;
      nFinished++;
    } else { // Need to look further
      snplist[pos] = rare;
    }
    nVariant++;

  }
  notice("Finished first pass. Processed %d rare variants across %d samples; finished %d.", nVariant, nsamples, nFinished);
  notice("%d pairs miss left ibs0; %d pairs miss right ibs0.", idpair_l.size(), idpair_r.size());
  if (nVariant < 1) {
    notice("Not enough rare variants in the given region");
    odw.close();
    return 0;
  }
  int32_t bound1 = (*(ibs0finder.posvec_que[0]))[0];
  int32_t bound2 = ibs0finder.posvec_que.back()->back();
  // Look backward
  int32_t wed = bound2;
  int32_t wst = bound1;
  int32_t pct = 0, fin = 0;
  std::string wreg;
  int32_t leftmost = pgmap.minpos;
  if (start > pgmap.centromere_st)
    leftmost = pgmap.centromere_ed;
  while (idpair_l.size() > 0 && (int32_t) idpair_l.size() > leftover/2 && wed > leftmost) {
    pct = 0;
    wed = wst - 1;
    wst = std::max(wed - bp_limit, 0);
    wreg = chrom + ":" + std::to_string(wst) + "-" + std::to_string(wed);
    if (!ibs0finder.Update_Fixed(wreg)) {
      continue;
    }
    auto itr = idpair_l.cbegin();
    while (itr != idpair_l.cend()) {
      // Iterate over pairs of individuals missing left ibs0 pt
      int32_t l = ibs0finder.FindIBS0((*itr).first.first,(*itr).first.second,wed,1);
      if (l > 0) {
        for (auto ptr : (*itr).second) {
          // Add this pos to all relevant variants
          fin = ptr->Add_half((*itr).first.first,(*itr).first.second,l,1);
          if (fin) {
            // Build tree
            ptr->Organize(pgmap);
            // Output this variant
            UpdateInfo_IBS0(odw.hdr, ptr->iv, ptr);
            odw.write(ptr->iv);
            // Delete the record
            snplist.erase(ptr->iv->pos+1);
            delete ptr;
            ptr = NULL;
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
  }
  if (wst <= leftmost) {
    auto itr = idpair_l.cbegin();
    while (itr != idpair_l.cend()) {
      int32_t l = leftmost;
      for (auto ptr : (*itr).second) {
        // Add this pos to all relevant variants
        fin = ptr->Add_half((*itr).first.first,(*itr).first.second,l,1);
        if (fin) {
          // Build tree
          ptr->Organize(pgmap);
          // Output this variant
          UpdateInfo_IBS0(odw.hdr, ptr->iv, ptr);
          odw.write(ptr->iv);
          // Delete the record
          snplist.erase(ptr->iv->pos+1);
          delete ptr;
          ptr = NULL;
          nFinished++;
        }
      }
      itr = idpair_l.erase(itr);
    }
  }

  // Look forward
  wed = bound2;
  wst = bound1;
  int32_t rightmost = pgmap.maxpos;
  if (end < pgmap.centromere_ed)
    rightmost = pgmap.centromere_st;
  while (idpair_r.size() > 0 && (int32_t) snplist.size() > leftover && wst < rightmost) {
    pct = 0;
    wst = wed + 1;
    wed = wst + bp_limit;
    wreg = chrom + ":" + std::to_string(wst) + "-" + std::to_string(wed);
    if (!ibs0finder.Update_Fixed(wreg)) {
      continue;
    }
    auto itr = idpair_r.cbegin();
    while (itr != idpair_r.cend()) {
      // Iterate over pairs of individuals missing left ibs0 pt
      int32_t r = ibs0finder.FindIBS0((*itr).first.first,(*itr).first.second,wst,0);
      if (r > 0) {
        for (auto ptr : (*itr).second) {
          // Add this pos to all relevant variants
          fin = ptr->Add_half((*itr).first.first,(*itr).first.second,r,2);
          if (fin) {
            // Build tree
            ptr->Organize(pgmap);
            // Output this variant
            UpdateInfo_IBS0(odw.hdr, ptr->iv, ptr);
            odw.write(ptr->iv);
            // Delete the record
            snplist.erase(ptr->iv->pos+1);
            delete ptr;
            ptr = NULL;
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
  }
  if (wst >= rightmost) {
    auto itr = idpair_r.cbegin();
    while (itr != idpair_r.cend()) {
      int32_t r = rightmost;
      for (auto ptr : (*itr).second) {
        // Add this pos to all relevant variants
        fin = ptr->Add_half((*itr).first.first,(*itr).first.second,r,2);
        if (fin) {
          // Build tree
          ptr->Organize(pgmap);
          // Output this variant
          UpdateInfo_IBS0(odw.hdr, ptr->iv, ptr);
          odw.write(ptr->iv);
          // Delete the record
          snplist.erase(ptr->iv->pos+1);
          delete ptr;
          ptr = NULL;
          nFinished++;
        }
      }
      itr = idpair_r.erase(itr);
    }
  }

  // Output those unfinished variants
  if (snplist.size() > 0) {
    for (auto ptr : snplist)
      odw.write(ptr.second->iv);
  }


  odw.close();

  return 0;
}









int32_t tmpUpdateCMInfo(int32_t argc, char** argv) {
  std::string inVcf, inMap, reg;
  std::string out;
  int32_t verbose = 100000;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("map",&inMap, "Map file for genetic distance")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")

  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  bp2cmMap pgmap(inMap, " ", "", -1, -1);
  notice("Read map. min %d; max %d; cst %d; ced %d.", pgmap.minpos, pgmap.maxpos, pgmap.centromere_st, pgmap.centromere_ed);

  std::vector<GenomeInterval> intervals;
  if (!reg.empty())
    parse_intervals(intervals, "", reg);
  BCFOrderedReader* odr = new BCFOrderedReader(inVcf, intervals);
  bcf1_t* iv = bcf_init();

  BCFOrderedWriter odw(out.c_str(),0);
  bcf_hdr_t* hnull = bcf_hdr_subset(odr->hdr, 0, 0, 0);
  bcf_hdr_remove(hnull, BCF_HL_FMT, NULL);
  odw.set_hdr(hnull);
  odw.write_hdr();

  int32_t *info_bp = NULL;
  float *info_cm = NULL;
  int32_t n_bp = 0, n_cm = 0;
  int32_t nVariants = 0;
  for(int32_t k = 0; odr->read(iv); ++k) {  // read marker
    if ( k % verbose == 0 )
      notice("Processing %d markers at %s:%d, updated %d variants.", k, bcf_hdr_id2name(odr->hdr, iv->rid), iv->pos+1, nVariants);

    bcf_unpack(iv, BCF_UN_ALL);
    bool update = 1;

    if (bcf_get_info_float(odr->hdr, iv, "AvgDist_cM", &info_cm, &n_cm) < 0) {update=0;}
    if (info_cm[0] > 0) {update=0;}
    if (bcf_get_info_int32(odr->hdr, iv, "AvgDist_bp", &info_bp, &n_bp) < 0) {update=0;}

    if (update) {

      int32_t st = iv->pos - info_bp[0]/2;
      int32_t ed = iv->pos + info_bp[0]/2;
      float newcm = pgmap.bpinterval2cm(st,ed);
      bcf1_t* nv = bcf_dup(iv);
      bcf_unpack(nv, BCF_UN_ALL);
      bcf_update_info_float(odw.hdr, nv, "AvgDist_cM", &newcm, 1);
      odw.write(nv);
      bcf_destroy(nv);
      nVariants++;
    } else {
      odw.write(iv);
    }

  }

  odw.close();
  return 0;
}

















