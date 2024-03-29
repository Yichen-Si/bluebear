#include "ibs0.h"

// Find shared rare allele for a pair of individual in one block
int32_t IndivRare::RareSharingFind(int32_t i, int32_t j,
                                   int32_t start, bool reverse) {

  if (rarelist[i].size() == 0 || rarelist[j].size() == 0) {return -1;}
  if (reverse) {
    auto pi = rarelist[i].rbegin();
    auto pj = rarelist[j].rbegin();
    while (pj != rarelist[j].rend() && (*pj) > start) {
      pj++;
    }
    while (pi != rarelist[i].rend() && (*pi) > start) {
      pi++;
    }
    while (pi != rarelist[i].rend() && pj != rarelist[j].rend()) {
      if ((*pj) > (*pi)) {
        pj++;
      } else if ((*pi) > (*pj)) {
        pi++;
      } else {
        return (*pi);
      }
    }
    return -1;
  } else {
    auto pi = rarelist[i].begin();
    auto pj = rarelist[j].begin();
    while (pj != rarelist[j].end() && (*pj) < start) {
      pj++;
    }
    while (pi != rarelist[i].end() && (*pi) < start) {
      pi++;
    }
    while (pi != rarelist[i].end() && pj != rarelist[j].end()) {
      if ((*pj) < (*pi)) {
        pj++;
      } else if ((*pi) < (*pj)) {
        pi++;
      } else {
        return (*pi);
      }
    }
    return -1;
  }
}


// Read one block, storing homozygote sites & rare allele carriers
int32_t ReadRareIBS0Block(std::string& reg, std::string& inVcf,
                   std::vector<bitmatrix*>& bmatRR_que,
                   std::vector<bitmatrix*>& bmatAA_que,
                   std::vector<IndivRare*>& rare_que,
                   std::vector<std::vector<int32_t>* >& posvec_que,
                   int32_t max_ac, int32_t min_hom_gts, bool add_to_front) {

  std::vector<GenomeInterval> ibs0interval;
  parse_intervals(ibs0interval, "", reg);
  BCFOrderedReader ibs0odr(inVcf, ibs0interval);
  int32_t nsamples = bcf_hdr_nsamples(ibs0odr.hdr);
  bcf1_t* ibs0iv = bcf_init();
  bitmatrix *bmatRR = new bitmatrix(nsamples);
  bitmatrix *bmatAA = new bitmatrix(nsamples);
  IndivRare *rmat   = new IndivRare(nsamples);
  std::vector<int32_t> *posvec = new std::vector<int32_t>;
  uint8_t *gtRR = (uint8_t*)calloc(nsamples, sizeof(uint8_t));
  uint8_t *gtAA = (uint8_t*)calloc(nsamples, sizeof(uint8_t));
  int32_t *p_gt = NULL;
  int32_t  n_gt = 0;
  int32_t *info_ac = NULL, *info_an = NULL;
  int32_t  n_ac = 0, n_an = 0;
  for(int32_t k=0; ibs0odr.read(ibs0iv); ++k) {  // read markers
    if ( bcf_get_genotypes(ibs0odr.hdr, ibs0iv, &p_gt, &n_gt) < 0 ) {
      notice("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(ibs0odr.hdr, ibs0iv->rid), ibs0iv->pos+1);
      return 0;
    }
    if (bcf_get_info_int32(ibs0odr.hdr, ibs0iv, "AC", &info_ac, &n_ac) < 0) {continue;}
    if (bcf_get_info_int32(ibs0odr.hdr, ibs0iv, "AN", &info_an, &n_an) < 0) {continue;}
    memset(gtRR, 0, nsamples);
    memset(gtAA, 0, nsamples);
    int32_t gcs[3] = {0,0,0};
    for(int32_t i=0; i < nsamples; ++i) {
      int32_t g1 = p_gt[2*i];
      int32_t g2 = p_gt[2*i+1];
      int32_t geno;
      if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
        //geno = 0;
      } else {
        geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
        if (info_ac[0] > info_an[0]/2) {
          geno = 2 - geno;
        }
        if ((info_ac[0] <= max_ac || info_an[0] - info_ac[0] <= max_ac) && geno > 0) {
          rmat->RareCarryAdd(i, ibs0iv->pos+1);
        }
        if ( geno == 0 )    { gtRR[i] = 1; }
        else if ( geno == 2 ) { gtAA[i] = 1; }
        ++gcs[geno];
      }
    }
    if (( gcs[0] >= min_hom_gts ) && ( gcs[2] >= min_hom_gts )) {
      bmatRR->add_row_bytes(gtRR);
      bmatAA->add_row_bytes(gtAA);
      posvec->push_back(ibs0iv->pos+1);
    }
  }
  free(gtRR);
  free(gtAA);
  bmatRR->transpose();
  bmatAA->transpose();
  if (bmatRR->ncol == 0) {
    delete bmatRR;
    delete bmatAA;
    delete posvec;
    delete rmat;
    return 0;
  }

  if (add_to_front) {
    bmatRR_que.insert(bmatRR_que.begin(), bmatRR);
    bmatAA_que.insert(bmatAA_que.begin(), bmatAA);
    posvec_que.insert(posvec_que.begin(), posvec);
    rare_que.insert(rare_que.begin(), rmat);
  } else {
    bmatRR_que.push_back(bmatRR);
    bmatAA_que.push_back(bmatAA);
    posvec_que.push_back(posvec);
    rare_que.push_back(rmat);
  }
  return bmatRR->ncol;
}



// Initialize ibs0 & rare variant sharing look up
RareIBS0lookup::RareIBS0lookup(const std::string &_inVcf, const std::string &_reg,
                       bp2cmMap &_pgmap, double _margin,
                       int32_t _ac, int32_t _ck, int32_t _mh) :
  inVcf(_inVcf), reg(_reg), pgmap(_pgmap), margin_cM(_margin), ac_cut(_ac),
  chunksize(_ck), min_hom_gts(_mh) {

  int32_t start, end, mid;
  int32_t ibs_st, ibs_ed;
  std::vector<std::string> v;
  std::string chrom;
  split(v, ":-", reg);
  if (!str2int32(v[1], start) || !str2int32(v[2], end)) {
    error("Invalid region for initialize ibs0 lookup blocks.");
  }
  chrom = v[0];
  mid = (start + end) / 2;

  ibs_st = mid - chunksize/2;
  ibs_ed = ibs_st + chunksize - 1;

  while ( ibs_st > pgmap.minpos && pgmap.bp2cm(start) - pgmap.bp2cm(ibs_st) < margin_cM ) {
    ibs_ed = ibs_st - 1;
    ibs_st -= chunksize;
  } // Find the starting point for adding ibs0 lookup blocks
  ibs_st = std::max(ibs_st, pgmap.minpos);

  // Moving forward
  while (ibs_st < pgmap.maxpos && pgmap.bp2cm(ibs_st)-pgmap.bp2cm(end) < margin_cM) {
    reg = chrom + ":" + std::to_string(ibs_st) + "-" + std::to_string(ibs_ed);
    if (ReadRareIBS0Block(reg, inVcf, bmatRR_que, bmatAA_que, rare_que, posvec_que, ac_cut, min_hom_gts, 0) > 0) {
      start_que.push_back((*posvec_que.back())[0]);
    }
    ibs_st = ibs_ed + 1;
    ibs_ed +=chunksize;
  } // Finish initialize lookup blocks
}

RareIBS0lookup::RareIBS0lookup(const std::string &_inVcf, const std::string &_reg,
                       bp2cmMap &_pgmap, int32_t _margin,
                       int32_t _ac, int32_t _ck, int32_t _mh) :
  inVcf(_inVcf), reg(_reg), pgmap(_pgmap), margin_bp(_margin), ac_cut(_ac),
  chunksize(_ck), min_hom_gts(_mh) {

  by_bp = 1;
  int32_t start, end, mid;
  int32_t ibs_st, ibs_ed;
  std::vector<std::string> v;
  std::string chrom;
  split(v, ":-", reg);
  if (!str2int32(v[1], start) || !str2int32(v[2], end)) {
    error("Invalid region for initialize ibs0 lookup blocks.");
  }
  chrom = v[0];
  mid = (start + end) / 2;

  // To jump over centromere
  int32_t off_left = 0, off_right = 0;
  if (pgmap.overlap_centro(start,end)) {
    off_left = std::max(start - pgmap.centromere_st, 0);
    off_right = std::max(end - pgmap.centromere_ed, 0);
  }

  // Locate the left most point to start
  ibs_st = mid - chunksize/2;
  ibs_ed = ibs_st + chunksize - 1;
  while ( ibs_st > pgmap.minpos && start - ibs_st < margin_bp + off_left ) {
    ibs_ed = ibs_st - 1;
    ibs_st -= chunksize;
  }
  ibs_st = std::max(ibs_st, pgmap.minpos);

  // Moving forward
  while (ibs_st < pgmap.maxpos && ibs_st-end < margin_bp + off_right) {
    reg = chrom + ":" + std::to_string(ibs_st) + "-" + std::to_string(ibs_ed);
    if (ReadRareIBS0Block(reg, inVcf, bmatRR_que, bmatAA_que, rare_que, posvec_que, ac_cut, min_hom_gts, 0) > 0) {
      start_que.push_back((*posvec_que.back())[0]);
    }
    ibs_st = ibs_ed + 1;
    ibs_ed +=chunksize;
  } // Finish initialize lookup blocks
}


// Find shared rare allele for a pair of individuals starting from pos
int32_t RareIBS0lookup::FindRare (int32_t i, int32_t j, int32_t pos, bool reverse) {

  if ( (*posvec_que[0])[0] > pos || posvec_que.back()->back() < pos ) {
    return -1;
  }
  int32_t k = (int32_t) start_que.size() - 1;
  while (start_que[k] > pos && k > 0) {k--;}
  int32_t rare = rare_que[k]->RareSharingFind(i,j,pos,reverse);
  if (reverse) {
    while (rare < 0 && k > 0) {
      k--;
      rare = rare_que[k]->RareSharingFind(i,j,pos,reverse);
    }
  } else {
    while(rare < 0 && k < (int32_t) start_que.size() - 1) {
      k++;
      rare = rare_que[k]->RareSharingFind(i,j,pos,reverse);
    }
  }
  return rare;
}

// Find IBS0 for a pair of individuals starting from pos
int32_t RareIBS0lookup::FindIBS0 (int32_t i, int32_t j, int32_t pos, bool reverse) {

  if ( (*posvec_que[0])[0] > pos || posvec_que.back()->back() < pos ) {
    return -1;
  }
  int32_t k = (int32_t) start_que.size() - 1;
  while (start_que[k] > pos && k > 0) {k--;}
  int32_t ibs0 = IBS0inOneBlock(bmatRR_que[k], bmatAA_que[k], posvec_que[k],
                                i, j, reverse, pos);
  if (reverse) {
    while (ibs0 < 0 && k > 0) {
      k--;
      ibs0 = IBS0inOneBlock(bmatRR_que[k], bmatAA_que[k], posvec_que[k],
                            i, j, reverse);
    }
  } else {
    while(ibs0 < 0 && k < (int32_t) start_que.size() - 1) {
      k++;
      ibs0 = IBS0inOneBlock(bmatRR_que[k], bmatAA_que[k], posvec_que[k],
                            i, j, reverse);
    }
  }
  return ibs0;
}


// Update to delete & add blocks to cover a new region
int32_t RareIBS0lookup::Update(const std::string &_reg) {

  reg = _reg;
  int32_t start, end;
  int32_t ibs_st, ibs_ed;
  std::vector<std::string> v;
  std::string chrom;
  split(v, ":-", reg);
  if (!str2int32(v[1], start) || !str2int32(v[2], end)) {
    error("Invalid region for update ibs0 lookup blocks.");
  }
  chrom = v[0];

  if (by_bp) {
    int32_t margin = margin_bp;
    // Delete unnecessary blocks
    while (start - posvec_que[0]->back() > margin &&
           posvec_que.size() > 0) {
      delete posvec_que[0]; posvec_que.erase(posvec_que.begin());
      delete bmatRR_que[0]; bmatRR_que.erase(bmatRR_que.begin());
      delete bmatAA_que[0]; bmatAA_que.erase(bmatAA_que.begin());
      delete rare_que[0]; rare_que.erase(rare_que.begin());
      start_que.erase(start_que.begin());
    }
    while (start_que.back() - end > margin &&
           posvec_que.size() > 0) {
      delete posvec_que.back(); posvec_que.pop_back();
      delete bmatRR_que.back(); bmatRR_que.pop_back();
      delete bmatAA_que.back(); bmatAA_que.pop_back();
      delete rare_que.back(); rare_que.pop_back();
      start_que.pop_back();
    }

    // Add new ibs0 lookup blocks to tail
    ibs_ed = posvec_que.back()->back();
    while (ibs_ed - end < margin) {
      ibs_st = ibs_ed + 1;
      ibs_ed = ibs_st + chunksize;
      if (ibs_st > pgmap.centromere_st && ibs_ed < pgmap.centromere_ed) {
        ibs_st = ibs_ed + 1;
        ibs_ed = ibs_st + chunksize;
        continue;
      }
      if (ibs_st > pgmap.maxpos)
        break;
      reg = chrom + ":" + std::to_string(ibs_st) + "-" + std::to_string(ibs_ed);
      if (ReadRareIBS0Block(reg, inVcf, bmatRR_que, bmatAA_que, rare_que, posvec_que, ac_cut, min_hom_gts, 0) > 0){
        start_que.push_back((*posvec_que.back())[0]);
      }
    }
    // Add new ibs0 lookup blocks to head
    ibs_st = start_que[0];
    while (start - ibs_st < margin) {
      ibs_ed = ibs_st - 1;
      ibs_st = ibs_ed - chunksize;
      if (ibs_st > pgmap.centromere_st && ibs_ed < pgmap.centromere_ed) {
        ibs_ed = ibs_st - 1;
        ibs_st = ibs_ed - chunksize;
        continue;
      }
      if (ibs_ed < pgmap.minpos)
        break;
      reg = chrom + ":" + std::to_string(ibs_st) + "-" + std::to_string(ibs_ed);
      if (ReadRareIBS0Block(reg, inVcf, bmatRR_que, bmatAA_que, rare_que, posvec_que, ac_cut, min_hom_gts, 1) > 0){
        start_que.push_back((*posvec_que[0])[0]);
      }
    }
  }

  else {
    double margin = margin_cM;
    // Delete unnecessary blocks
    while (pgmap.bp2cm(start) - pgmap.bp2cm(posvec_que[0]->back()) > margin &&
           posvec_que.size() > 0) {
      delete posvec_que[0]; posvec_que.erase(posvec_que.begin());
      delete bmatRR_que[0]; bmatRR_que.erase(bmatRR_que.begin());
      delete bmatAA_que[0]; bmatAA_que.erase(bmatAA_que.begin());
      delete rare_que[0]; rare_que.erase(rare_que.begin());
      start_que.erase(start_que.begin());
    }
    while (pgmap.bp2cm((*posvec_que.back())[0]) - pgmap.bp2cm(end) > margin &&
           posvec_que.size() > 0) {
      delete posvec_que.back(); posvec_que.pop_back();
      delete bmatRR_que.back(); bmatRR_que.pop_back();
      delete bmatAA_que.back(); bmatAA_que.pop_back();
      delete rare_que.back(); rare_que.pop_back();
      start_que.pop_back();
    }

    // Add new ibs0 lookup blocks to tail
    ibs_ed = posvec_que.back()->back();
    while (pgmap.bp2cm(ibs_ed) - pgmap.bp2cm(end) < margin) {
      ibs_st = ibs_ed + 1;
      ibs_ed = ibs_st + chunksize;
      if (ibs_st > pgmap.centromere_st && ibs_ed < pgmap.centromere_ed) {
        ibs_st = ibs_ed + 1;
        ibs_ed = ibs_st + chunksize;
        continue;
      }
      if (ibs_st > pgmap.maxpos)
        break;
      reg = chrom + ":" + std::to_string(ibs_st) + "-" + std::to_string(ibs_ed);
      if (ReadRareIBS0Block(reg, inVcf, bmatRR_que, bmatAA_que, rare_que, posvec_que, ac_cut, min_hom_gts, 0) > 0){
        start_que.push_back((*posvec_que.back())[0]);
      }
    }
    // Add new ibs0 lookup blocks to head
    ibs_st = (*posvec_que[0])[0];
    while (pgmap.bp2cm(start) - pgmap.bp2cm(ibs_st) < margin) {
      ibs_ed = ibs_st - 1;
      ibs_st = ibs_ed - chunksize;
      if (ibs_st > pgmap.centromere_st && ibs_ed < pgmap.centromere_ed) {
        ibs_ed = ibs_st - 1;
        ibs_st = ibs_ed - chunksize;
        continue;
      }
      if (ibs_ed < pgmap.minpos)
        break;
      reg = chrom + ":" + std::to_string(ibs_st) + "-" + std::to_string(ibs_ed);
      if (ReadRareIBS0Block(reg, inVcf, bmatRR_que, bmatAA_que, rare_que, posvec_que, ac_cut, min_hom_gts, 1) > 0){
        start_que.push_back((*posvec_que[0])[0]);
      }
    }
  }

  return (int32_t) start_que.size();

}

