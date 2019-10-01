#include "cramore.h"
#include "utils.h"
#include "bcf_filter_arg.h"
#include "bcf_ordered_reader.h"
#include "compact_matrix.h"
#include <iomanip>
#include "bp2cm.h"
#include "ibs0.h"

class IBS0rec {
public:
  int32_t st = 0, ed = 0;
  std::vector<int32_t> rarepos;
  std::vector<int32_t> rareac;

bool update(int32_t pos, bool front = 0) {
  if (front) {
    if (ed != 0) {
      return 1; // Complete
    } else {
      st = pos; // Missing ed pt
      return 0;
    }
  } else {
    if (rarepos.size() > 0 && st != 0) {
      return 1; // Complete
    } else if (rarepos.size() > 0 && ed == 0) {
      ed = pos; // Missing st pt
      return 0;
    } else if (rarepos.size() == 0) {
      st = pos; // Not start yet
      return 0;
    } else {
      return 0;
    }
  }
}

bool update(int32_t pos, int32_t ac) {
  rarepos.push_back(pos);
  rareac.push_back(ac);
  if (st != 0) {
    return 0;
  } else {
    return 1;
  }
}

void clear() {
  st = 0;
  rarepos.clear();
  rareac.clear();
}

};


struct interval {
  int32_t id1, id2;
  int32_t st, ed;
  double cm;
  interval(int32_t _x, int32_t _y, int32_t _st, int32_t _ed, double _l) :
  id1(_x), id2(_y), st(_st), ed(_ed), cm(_l) {}
};

class RareVariant {
public:
  int32_t ac;
  std::vector<interval> ibs0pair;
  int32_t ovst = -1, oved = -1; // st & ed of the overlap region

RareVariant(int32_t _ac) : ac(_ac) {}
bool Add(int32_t id1, int32_t id2, int32_t st, int32_t ed, double cm) {
  interval rec(id1,id2,st,ed,cm);
  ibs0pair.push_back(rec);
  if (ovst < 0 || oved < 0) {
    ovst = st; oved = ed;
  } else {
    ovst = std::max(ovst, st);
    oved = std::min(oved, ed);
  }
  if ( ((int32_t) ibs0pair.size()) >= (ac*(ac-1)/2) ) {
    return 1;
  } else {
    return 0;
  }
}
void Organize() {
  // Sort individual pairs by id, so the output can be re-shaped to a matrix
  std::sort(ibs0pair.begin(), ibs0pair.end(),
           [](const interval& a, const interval& b) {
           return a.id1 < b.id1;});
}

};

// goal -- get all pairwise no-ibs0 regions starting from the given region
//         covering shared rare variants
int32_t IBS0PairwiseScan(int32_t argc, char** argv) {
  std::string inVcf, inMap, chrom, reg, oreg;
  std::string out, outf;
  int32_t min_hom_gts = 1;
  int32_t verbose = 10000;
  int32_t detailed = 0;
  int32_t batch_size = 10000;
  int32_t min_variant = 1;
  int32_t max_rare_ac = 10;
  int32_t n_samples = 0;
  int32_t start, end;
  int32_t leftover = 20, maxoutreach = 3000000;
  int32_t cst = -1, ced = -1;

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
    LONG_INT_PARAM("batch-size",&batch_size, "Size of batches (in # of samples) to calculate the no-IBS0 pairs")
    LONG_INT_PARAM("min-variant",&min_variant, "Minimum number of variants to present to have output file")
    LONG_INT_PARAM("left-over",&leftover, "Stop when only x pairs are left")
    LONG_INT_PARAM("out-reach",&maxoutreach, "How far to look forward & backward")
    // For testing
    LONG_INT_PARAM("num-samples",&n_samples, "Number of samples to test")

    LONG_INT_PARAM("max-rare",&max_rare_ac, "Maximal minor allele count to be considered as anchor for IBD")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")
    LONG_INT_PARAM("detailed",&detailed,"If output all shared rare variatns POS:AC")


  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // sanity check of input arguments
  if ( inVcf.empty() || reg.empty() || inMap.empty() || out.empty() ) {
    error("[E:%s:%d %s] --in-vcf, --region, --inMap, --out are required parameters",__FILE__,__LINE__,__FUNCTION__);
  }

  // Genetic map
  bp2cmMap pgmap(inMap, " ", cst, ced);

  // Output
  std::ofstream wf;
  outf = out + "_" + std::to_string(max_rare_ac) + "_" + reg + ".list";
  wf.open(outf, std::ios::trunc);

  // Region to process
  oreg = reg;
  std::vector<std::string> v;
  split(v, ":-", reg);
  chrom = v[0];
  if (!str2int32(v[1], start) || !str2int32(v[2], end)) {
    error("Invalid region.");
  }

  // bcf reader
  std::vector<GenomeInterval> intervals;
  parse_intervals(intervals, "", reg);
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();

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

  int32_t nVariant = 0, nMissST = 0;
  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);

  notice("Started Reading site information from VCF file, identifying %d samples", nsamples);
  if ( nsamples == 0 )
    error("FATAL ERROR: The VCF does not have any samples with genotypes");

  // For testing
  if (n_samples > 0 && nsamples > n_samples)
    nsamples = n_samples;

  bitmatrix bmatRR(nsamples);
  bitmatrix bmatAA(nsamples);
  uint8_t* gtRR = (uint8_t*)calloc(nsamples, sizeof(uint8_t));
  uint8_t* gtAA = (uint8_t*)calloc(nsamples, sizeof(uint8_t));
  std::vector<int32_t> mac; // minor allele count
  std::vector<int32_t> snoopy; // informative position

  std::map<int32_t, std::map<int32_t, IBS0rec*>> idpair;
  std::map<int32_t, RareVariant*> snplist;
  int32_t* p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t *info_ac = NULL;
  int32_t *info_an = NULL;
  int32_t n_ac = 0, n_an = 0;

  // First pass.
  // Read markers with both homozygotes. Find indiv. pairs sharing rare allele.
  for(int32_t k=0; odr.read(iv); ++k) {
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
    int32_t ifflip = 0, ac = 0, an = 0;
    ac = info_ac[0]; an = info_an[0];
    if (ac > an/2) {
      ac = an - ac;
      ifflip = 1;
    }
    if (ac < 2) {continue;}
    memset(gtRR, 0, nsamples);
    memset(gtAA, 0, nsamples);
    int32_t gcs[3] = {0,0,0};
    std::vector<int32_t> carriers;
    for(int32_t i=0; i < nsamples; ++i) {
      int32_t g1 = p_gt[2*i];
      int32_t g2 = p_gt[2*i+1];
      int32_t geno;
      if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
         //geno = 0;
      }
      else {
        geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
        if ( ifflip == 1 ) {geno = 2 - geno;}
        if ( geno == 0 )    { gtRR[i] = 1; }
        else if ( geno == 2 ) { gtAA[i] = 1; }
        ++gcs[geno];
        if ( ac <= max_rare_ac && geno == 1 ) {
          // If it is a rare variant, 1 labels carriers
          gtAA[i] = 1; gtRR[i] = 1;
          carriers.push_back(i);
        }
      }
    }
    if ( (ac > max_rare_ac) && (( gcs[0] < min_hom_gts ) || ( gcs[2] < min_hom_gts )) ) {
      continue;
    } else if (ac <= max_rare_ac && ac != gcs[1]+gcs[2]) {
      continue;
    } else {
      bmatRR.add_row_bytes(gtRR);
      bmatAA.add_row_bytes(gtAA);
      mac.push_back(ac);
      snoopy.push_back(iv->pos+1);
      if (ac <= max_rare_ac) {
        ++nVariant;
        snplist[iv->pos+1] = new RareVariant(ac);
        for (uint32_t i = 0; i < carriers.size()-1; ++i) {
          for (uint32_t j = i+1; j < carriers.size(); ++j) {
            idpair[carriers[i]][carriers[j]] = NULL;
          }
        }
      }
    }
  }
  notice("Finished first pass. Read %d SNPs, including %d rare variants across %d samples.", mac.size(), nVariant, nsamples);

  free(gtRR);
  free(gtAA);

  if ( nVariant < min_variant ) {
    notice("Observed only %d rare variants. Skipping IBD segment detection for this chunk...", nVariant);
    return 0;
  }

  bmatRR.transpose();
  bmatAA.transpose();

  // Initialize pairwise recording
  int32_t pairct = 0;
  for (auto & k : idpair) {
    for (auto & k2 : k.second) {
      k2.second = new IBS0rec;
    }
    pairct += k.second.size();
  }

  notice("Searching for IBS0 segments among %d pairs.", pairct);

  for (auto & k1 : idpair) {
    int32_t i = k1.first;
    uint8_t* iRR = bmatRR.get_row_bits(i);
    uint8_t* iAA = bmatAA.get_row_bits(i);
    for (auto & k2 : k1.second) {
      int32_t j = k2.first;
      uint8_t* jRR = bmatRR.get_row_bits(j);
      uint8_t* jAA = bmatAA.get_row_bits(j);
      IBS0rec * v = k2.second; // Only for clearity
      for(int32_t k=0; k < bmatRR.nbytes_col; ++k) {
        uint8_t byte = ( iRR[k] ^ jRR[k] ) & ( iAA[k] ^ jAA[k] );
        uint8_t rarebyte = ( iRR[k] & jRR[k] ) & ( iAA[k] & jAA[k] );
        if (!byte && !rarebyte) {continue;}
        for (int32_t bit=0; bit<8 && k*8+bit<bmatRR.ncol; ++bit) {
          int32_t pt = k*8+bit;
          if (mac[pt] <= max_rare_ac && ((rarebyte >> (7-bit)) & 0x01) ) {
          // share rare allele
            if (v->update(snoopy[pt], mac[pt])) {nMissST++;}
          }
          else if ( (byte >> (7-bit)) & 0x01 ) {
          // IBS0
            if (v->update(snoopy[pt])) {
            // End of an informative no-ibs0 region
              double lcm = pgmap.bp2cm(snoopy[pt])-pgmap.bp2cm(v->st);
              for (uint32_t it = 0; it < (v->rarepos).size(); ++it) {
                const auto & mtr = snplist.find(v->rarepos[it]);
                if (mtr != snplist.end()) {
                  RareVariant * ptr = mtr->second;
                  if (ptr->Add(i,j,v->st,snoopy[pt],lcm)) {
                    // Output this rare variant
                    ptr->Organize();
                    wf << v->rarepos[it] << '\t' << v->rareac[it] << '\t';
                    wf << ptr->ovst << '\t' << ptr->oved << '\t';
                    wf << std::fixed << std::setprecision(3) <<
                    pgmap.bp2cm(ptr->oved)-pgmap.bp2cm(ptr->ovst) << '\t';
                    for (auto & x : ptr->ibs0pair) {
                      wf << std::fixed << std::setprecision(3) << x.cm << ',';
                    }
                    wf << '\n';
                    delete ptr;
                    snplist.erase(v->rarepos[it]);
                  }
                }
              }
              v->clear();
            }
          }
        }
      }
    } // Finish updating one pair of individuals
  } // Finish updating all pairs (that share >= 1 rare alleles)

  notice("Finish processing within the given region. %d variants left.", snplist.size());
  if (bmatRR.bytes != NULL) {free(bmatRR.bytes); bmatRR.bytes=NULL;}
  if (bmatAA.bytes != NULL) {free(bmatAA.bytes); bmatAA.bytes=NULL;}
  mac.clear();

  pairct = 0;
  int32_t missct = 0;
  for (auto k1 = idpair.begin(); k1 != idpair.end(); ) {
    for (auto k2 = k1->second.begin(); k2 != k1->second.end(); ) {
      // if (k2->second->st==0 && (k2->second->rarepos).size() == 0) {
      if ((k2->second->rarepos).size() == 0) {
        delete k2->second;
        k1->second.erase(k2++);
      } else {
        if (k2->second->st == 0) {missct++;}
        k2++;
      }
    }
    if (k1->second.size() == 0) {
      idpair.erase(k1++);
    } else {
      pairct += k1->second.size();
      k1++;
    }
  } // Delete finished no-ibs0 records.

  notice("%d no-ibs0 regions remain open; %d(%d) regions miss start points.", pairct,missct,nMissST);

  // Look backward. Find starting points of existing no-ibs0 regions.
  int32_t win = 200000;
  int32_t wed = start-1;
  int32_t wst = wed - win + 1;
  double maf = 0.005;
  while (wst > pgmap.minpos && missct > leftover && start - wst < maxoutreach) {
    reg = chrom + ":" + std::to_string(wst) + "-" + std::to_string(wed);
    parse_intervals(intervals, "", reg);
    odr.jump_to_interval(intervals[0]);
    bcf_clear(iv);
    bitmatrix bmatRR(nsamples);
    bitmatrix bmatAA(nsamples);
    snoopy.clear(); // informative position
    uint8_t* gtRR = (uint8_t*)calloc(nsamples, sizeof(uint8_t));
    uint8_t* gtAA = (uint8_t*)calloc(nsamples, sizeof(uint8_t));
    // Read markers with both homozygotes.
    for(int32_t k=0; odr.read(iv); ++k) {
      if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
      if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
      double af = info_ac[0]*1.0/info_an[0];
      if (af < maf || af > 1-maf) {continue;}
      // extract genotype and apply genotype level filter
      if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
        error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
      }
      memset(gtRR, 0, nsamples);
      memset(gtAA, 0, nsamples);
      int32_t gcs[3] = {0,0,0};
      for(int32_t i=0; i < nsamples; ++i) {
        int32_t g1 = p_gt[2*i];
        int32_t g2 = p_gt[2*i+1];
        int32_t geno = 0;
        if ( !bcf_gt_is_missing(g1) && !bcf_gt_is_missing(g2) ) {
          geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
          if ( geno == 0 )    { gtRR[i] = 1; }
          else if ( geno == 2 ) { gtAA[i] = 1; }
          ++gcs[geno];
        }
      }
      if (( gcs[0] < min_hom_gts ) || ( gcs[2] < min_hom_gts )) {
        continue;
      } else {
        bmatRR.add_row_bytes(gtRR);
        bmatAA.add_row_bytes(gtAA);
        snoopy.push_back(iv->pos+1);
      }
    }
    free(gtRR);
    free(gtAA);
    notice("Second stage. Processing %s, Read %d SNPs. Still missing %d.", reg.c_str(), snoopy.size(), missct);
    if ( (int32_t) snoopy.size() < 10 ) {
      notice("Observed only %d variants. Skipping IBD segment detection for this chunk...", snoopy.size());
      wed = wst-1;
      wst = wed-win+1;
      continue;
    }

    bmatRR.transpose();
    bmatAA.transpose();
    for (auto & k1 : idpair) {
      int32_t i = k1.first;
      for (auto & k2 : k1.second) {
        if (k2.second->st != 0) {continue;}
        int32_t j = k2.first;
        int32_t previbs0 = IBS0inOneBlock(&bmatRR,&bmatAA,&snoopy,i,j,1);
        if (previbs0 < 0) {continue;}
        missct--;
        if (k2.second->update(previbs0, (bool) 1)) {
          pairct--;
        // Completed an informative no-ibs0 region
          IBS0rec * v = k2.second; // Only for clearity
          double lcm = pgmap.bp2cm(v->ed)-pgmap.bp2cm(previbs0);
          for (uint32_t it = 0; it < (v->rarepos).size(); ++it) {
            const auto & mtr = snplist.find(v->rarepos[it]);
            if (mtr != snplist.end()) {
              RareVariant * ptr = mtr->second;
              // RareVariant * ptr = snplist[v->rarepos[it]];
              if (ptr->Add(i,j,previbs0,v->ed,lcm)) {
                // Output this rare variant
                ptr->Organize();
                wf << v->rarepos[it] << '\t' << v->rareac[it] << '\t';
                wf << ptr->ovst << '\t' << ptr->oved << '\t';
                wf << std::fixed << std::setprecision(3) <<
                pgmap.bp2cm(ptr->oved)-pgmap.bp2cm(ptr->ovst) << '\t';
                for (auto & x : ptr->ibs0pair) {
                  wf << std::fixed << std::setprecision(3) << x.cm << ',';
                }
                wf << '\n';
                delete ptr;
                snplist.erase(v->rarepos[it]);
              }
            }
          }
          delete k2.second;
          k1.second.erase(k2.first);
        }
      } // Finish updating one pair of individuals
    } // Finish updating all pairs (that share >= 1 rare alleles)
    wed = wst-1;
    wst = wed-win+1;
  }

  pairct = 0;
  missct = 0;
  for (auto k1 = idpair.begin(); k1 != idpair.end(); ) {
    for (auto k2 = k1->second.begin(); k2 != k1->second.end(); ) {
      // if (k2->second->st==0 && (k2->second->rarepos).size() == 0) {
      if ((k2->second->rarepos).size() == 0) {
        delete k2->second;
        k1->second.erase(k2++);
      } else {
        if (k2->second->st == 0) {missct++;}
        k2++;
      }
    }
    if (k1->second.size() == 0) {
      idpair.erase(k1++);
    } else {
      pairct += k1->second.size();
      k1++;
    }
  } // Delete finished no-ibs0 records.
  notice("%d no-ibs0 regions remain open; %d regions miss start points.", pairct,missct);



  // Finding the ends of the existing no-ibs0 regions
  pairct -= missct;
  reg = chrom + ":" + std::to_string(end) + "-" + std::to_string(pgmap.maxpos);
  parse_intervals(intervals, "", reg);
  odr.jump_to_interval(intervals[0]);
  bcf_clear(iv);

  bool flag = 0;
  for(int32_t k=0; odr.read(iv); ++k) {
    if (iv->pos - end > maxoutreach) {
      break;
    }
    if ( k % verbose == 0 )
      notice("Final stage: processing %d markers at %s:%d; %d (%d) pairs left.", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, pairct, idpair.size());
    if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
    }
    if (bcf_get_info_int32(odr.hdr, iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (bcf_get_info_int32(odr.hdr, iv, "AN", &info_an, &n_an) < 0) {continue;}
    double af = info_ac[0]*1.0/info_an[0];
    if (af < maf || af > 1-maf) {continue;}
    for (auto k1 = idpair.begin(); k1 != idpair.end(); ) {
      for (auto k2 = k1->second.begin(); k2 != k1->second.end(); ) {
        if (k2->second->st == 0) {
          k2++;
          continue;
        }
        int32_t id1 = k1->first, id2 = k2->first;
        if ( bcf_gt_is_missing(p_gt[2*id1]) || bcf_gt_is_missing(p_gt[2*id1+1]) || bcf_gt_is_missing(p_gt[2*id2]) || bcf_gt_is_missing(p_gt[2*id2+1]) ) {
          k2++; continue;
        }
        int32_t g1 = ((bcf_gt_allele(p_gt[2*id1]) > 0) ? 1 : 0) + ((bcf_gt_allele(p_gt[2*id1+1]) > 0) ? 1 : 0);
        int32_t g2 = ((bcf_gt_allele(p_gt[2*id2]) > 0) ? 1 : 0) + ((bcf_gt_allele(p_gt[2*id2+1]) > 0) ? 1 : 0);
        if (g1*g2==0 && g1+g2==2) {
          if (k2->second->update(iv->pos+1)) {
            // End of an informative no-ibs0 region
            IBS0rec * v = k2->second;
            double lcm = pgmap.bp2cm(iv->pos+1)-pgmap.bp2cm(v->st);
            for (uint32_t it = 0; it < (v->rarepos).size(); ++it) {
              const auto & mtr = snplist.find(v->rarepos[it]);
              if (mtr != snplist.end()) {
                RareVariant * ptr = mtr->second;
                if (ptr->Add(id1,id2,v->st,iv->pos+1,lcm)) {
                  // Output this rare variant
                  ptr->Organize();
                  wf << v->rarepos[it] << '\t' << v->rareac[it] << '\t';
                  wf << ptr->ovst << '\t' << ptr->oved << '\t';
                  wf << std::fixed << std::setprecision(3) <<
                  pgmap.bp2cm(ptr->oved)-pgmap.bp2cm(ptr->ovst) << '\t';
                  for (auto & x : ptr->ibs0pair) {
                    wf << std::fixed << std::setprecision(3) << x.cm << ',';
                  }
                  wf << '\n';
                  delete ptr;
                  snplist.erase(v->rarepos[it]);
                }
              }
            }
            delete k2->second;
            k1->second.erase(k2++);
            pairct--;
          } else {k2++;}
        } else {k2++;}
      }
      if (k1->second.size() == 0) {
        idpair.erase(k1++);
        if (pairct < leftover) {
          flag = 1; break;
        }
      } else {k1++;}
      if (flag) {break;}
    }
    if (flag) {break;}
  }

  notice("Finished searching for potential IBD segments covering shared rare allele with AC <= %d", max_rare_ac);
  wf.close();

  notice("%d variants left.", snplist.size());

  outf = out + "_" + std::to_string(max_rare_ac) + "_" + oreg + ".left.pair";
  wf.open(outf, std::ios::trunc);
  for (auto & k1 : idpair) {
    if (k1.second.size() <= 0) {continue;}
    for (auto & k2 : k1.second) {
      wf << k1.first << '\t' << k2.first << '\t' << k2.second->st << '\t' << k2.second->ed << '\t';
      for (auto & v : k2.second->rarepos) {
        wf << v << ',';
      }
      wf << '\n';
    }
  }
  wf.close();

  outf = out + "_" + std::to_string(max_rare_ac) + "_" + oreg + ".left.snp";
  wf.open(outf, std::ios::trunc);
  for (auto & v : snplist) {
    RareVariant * ptr = v.second;
    ptr->Organize();
    wf << v.first << '\t' << ptr->ac << '\t';
    wf << ptr->ovst << '\t' << ptr->oved << '\t';
    for (auto & x : ptr->ibs0pair)
      wf << x.id1 << ',' << x.id2 << ',' << std::fixed << std::setprecision(3) << x.cm << ';';
    wf << '\n';
  }
  wf.close();

  return 0;
}

