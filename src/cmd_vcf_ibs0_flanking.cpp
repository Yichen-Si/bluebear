#include "cramore.h"
#include "utils.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "ibs0.h"
#include "rare_variant_config.h"
#include <iomanip>

class ibs0pair {
  public:
    int32_t pos, ac;
    std::string id1, id2;
    int32_t left_bp=-1, right_bp=-1;
    double length_cM;
    std::string info;

  ibs0pair(int32_t _p, int32_t _ac, std::string _v=".") : pos(_p), ac(_ac), info(_v) {}
  void Add_id(std::string i, std::string j) {
    id1 = i;
    id2 = j;
  }
  void Add_info(std::string& v) {info = v;}
  void Add_half(int32_t pt, int32_t direction) {
    if (direction == 1) {
      left_bp = pt;
    } else {
      right_bp = pt;
    }
  }
  bool Finished() {
    return (left_bp > 0 && right_bp > 0);
  }
};

void output_pair(FILE * wf, bp2cmMap& pgmap, ibs0pair* snp, std::string& chrom) {

  double cm = pgmap.bpinterval2cm( snp->left_bp, snp->right_bp );
  fprintf(wf, "%s\t%d\t%d\t%s\t%s\t%d\t%d\t%f\t%s\n",
          chrom.c_str(), snp->pos, snp->ac, snp->id1.c_str(), snp->id2.c_str(),
          snp->left_bp, snp->right_bp, cm, snp->info.c_str() );
}

// goal -- Input option 1: a list of pos & individual pairs
//         Input option 2: a AC criterion
//         Output: flanking IBS0 information for individual pairs

int32_t cmdVcfIBS0Flank(int32_t argc, char** argv) {
  std::string inVcf, inQuery, inMap, chrom, reg, oreg;
  std::string out;
  int32_t min_hom_gts = 1;
  int32_t verbose = 1000;
  int32_t min_variant = 1;
  int32_t min_rare_ac = 2;
  int32_t max_rare_ac = 2;
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
    LONG_STRING_PARAM("in-query",&inQuery, "Position and individual pairs to query")

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
    LONG_INT_PARAM("min-ac",&min_rare_ac, "Minimum minor allele count to be considered as anchor for IBD")
    LONG_INT_PARAM("max-ac",&max_rare_ac, "Maximum minor allele count to be considered as anchor for IBD")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output file name")
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

  // output
  FILE * wf;
  wf = fopen(out.c_str(), "w");
  fputs("CHR\tPOS\tAC\tID1\tID2\tLeft\tRight\tLength_cM\tInfo\n", wf);

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
  char** id_ptr = bcf_hdr_get_samples(odr.hdr);
  std::vector<std::string> id_samples(id_ptr, id_ptr+nsamples);
  std::map<std::string, int32_t> id_index;
  for (int32_t i = 0; i < nsamples; ++i)
    id_index[id_samples[i]] = i;

  notice("Started Reading site information from VCF file, identifying %d samples", nsamples);
  if ( nsamples == 0 )
    error("FATAL ERROR: The VCF does not have any samples with genotypes");

  IBS0lookup ibs0finder(inVcf, reg, pgmap, bp_limit/2, ck_len, 1);
  if (ibs0finder.start_que.size() < 1) {
    notice("Not enough variants in the given region. Stopped without output.");
    return 0;
  }

  std::map< std::pair<int32_t, int32_t>, std::vector<ibs0pair*> > idpair_l;
  std::map< std::pair<int32_t, int32_t>, std::vector<ibs0pair*> > idpair_r;
  std::vector<ibs0pair*> snplist;
  int32_t* p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t *info_ac = NULL;
  int32_t *info_an = NULL;
  int32_t n_ac = 0, n_an = 0, pos = 0;


    // First pass - get pairs & focal pos to look for IBS0
  if (inQuery.empty()) {

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
      if (ac < min_rare_ac || ac > max_rare_ac) {continue;}

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

      // Creat a new ibs0pair record
      pos = iv->pos+1;
        // Find pairwise ibs0 w/in a short limit
      for (int32_t i = 0; i < ac-1; ++i) {
        for (int32_t j = i+1; j < ac; ++j) {
          ibs0pair* rare = new ibs0pair(pos, ac);
          rare->Add_id(id_samples[carry[i]], id_samples[carry[j]]);
          int32_t r = ibs0finder.FindIBS0(carry[i],carry[j],pos,0);
          int32_t l = ibs0finder.FindIBS0(carry[i],carry[j],pos,1);
          if (l > 0 && r > 0) {
            rare->Add_half(l, 1);
            rare->Add_half(r, 0);
            output_pair(wf, pgmap, rare, chrom);
            delete rare;
            nFinished++;
          } else {
            std::pair<int32_t, int32_t> kpair(carry[i],carry[j]);
            if (l > 0) {
              rare->Add_half(l, 1);
              idpair_r[kpair].push_back(rare);
            } else if (r > 0) {
              rare->Add_half(r, 0);
              idpair_l[kpair].push_back(rare);
            } else {
              idpair_r[kpair].push_back(rare);
              idpair_l[kpair].push_back(rare);
            }
            snplist.push_back(rare);
          }
        }
      }
      nVariant++;
    }

  } else {

    std::ifstream ifs;
    std::string line;
    std::vector<std::string> words;
    std::vector<int32_t> idvec(2,0);
    ifs.open(inQuery, std::ifstream::in);
    if (!ifs.is_open()) {
      error("Query file cannot be opened");
    }
    while(std::getline(ifs, line)) {
      words.clear();
      split(words, "\t", line);
      auto ptr1 = id_index.find(words[3]);
      auto ptr2 = id_index.find(words[4]);
      if (ptr1==id_index.end() || ptr2==id_index.end()) {continue;}
      int32_t p, a;
      if (!str2int32(words[1], p)) {continue;}
      if (!str2int32(words[2], a)) {continue;}
      if (p < start) {continue;}
      if (p >= end) {break;}
      ibs0pair* rare = new ibs0pair(p, a);
      rare->Add_id(words[3], words[4]);
      if (words.size() > 5) {
        rare->Add_info(words[5]);
      }
      int32_t r = ibs0finder.FindIBS0(ptr1->second,ptr2->second,p,0);
      int32_t l = ibs0finder.FindIBS0(ptr1->second,ptr2->second,p,1);
      if (l > 0 && r > 0) {
        rare->Add_half(l, 1);
        rare->Add_half(r, 0);
        output_pair(wf, pgmap, rare, chrom);
        delete rare;
        nFinished++;
      } else {
        std::pair<int32_t, int32_t> kpair(ptr1->second,ptr2->second);
        if (l > 0) {
          rare->Add_half(l, 1);
          idpair_r[kpair].push_back(rare);
        } else if (r > 0) {
          rare->Add_half(r, 0);
          idpair_l[kpair].push_back(rare);
        } else {
          idpair_r[kpair].push_back(rare);
          idpair_l[kpair].push_back(rare);
        }
        snplist.push_back(rare);
      }
    }
    nVariant = snplist.size();
    ifs.close();
  }

  notice("Finished first pass. Processed %d rare variants across %d samples; finished %d.", nVariant, nsamples, nFinished);
  notice("%d pairs miss left ibs0; %d pairs miss right ibs0.", idpair_l.size(), idpair_r.size());
  if (nVariant < 1) {
    notice("Not enough rare variants in the given region");
    fclose(wf);
    return 0;
  }
  int32_t bound1 = (*(ibs0finder.posvec_que[0]))[0];
  int32_t bound2 = ibs0finder.posvec_que.back()->back();
  nFinished = 0;

  // Look backward
  int32_t wed = bound2;
  int32_t wst = bound1;
  int32_t pct = 0;
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
        for (auto & ptr : (*itr).second) {
          // Add this pos to all relevant variants
          ptr->Add_half(l,1);
          if ( ptr->Finished() ) {
            output_pair(wf, pgmap, ptr, chrom);
            delete ptr;
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
      for (auto & ptr : (*itr).second) {
        // Add this pos to all relevant variants
        ptr->Add_half(l,1);
        if ( ptr->Finished() ) {
          output_pair(wf, pgmap, ptr, chrom);
          delete ptr;
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
        for (auto & ptr : (*itr).second) {
          // Add this pos to all relevant variants
          ptr->Add_half(r,2);
          if ( ptr->Finished() ) {
            output_pair(wf, pgmap, ptr, chrom);
            delete ptr;
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
      for (auto & ptr : (*itr).second) {
        // Add this pos to all relevant variants
        ptr->Add_half(r,2);
        if ( ptr->Finished() ) {
          output_pair(wf, pgmap, ptr, chrom);
          delete ptr;
          nFinished++;
        }
      }
      itr = idpair_r.erase(itr);
    }
  }

  fclose(wf);

  return 0;
}







