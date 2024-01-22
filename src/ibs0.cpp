#include "ibs0.h"

// Find the first IBS0 position for a pair of individuals
// Return -1 if no IBS0 is found, -2 if error
int32_t IBS0inOneBlock(bitmatrix* bmatRR, bitmatrix* bmatAA,
                       std::vector<int32_t>* posvec,
                       int32_t i, int32_t j, bool reverse,
                       int32_t start, int32_t _ptr) {

  uint8_t* iRR = bmatRR->get_row_bits(i);
  uint8_t* iAA = bmatAA->get_row_bits(i);
  uint8_t* jRR = bmatRR->get_row_bits(j);
  uint8_t* jAA = bmatAA->get_row_bits(j);
  int32_t  ptr = 0; // Index in this block
  int32_t  k, bpos; // Index in bitmat

  // Locate the first position <= start
  if (_ptr >= 0) {
    ptr = _ptr;
  } else {
    ptr = posvec->size() - 1;
    if (start > 0) {
      while(ptr > 0 && (*posvec)[ptr] > start) {ptr--;}
    } else {
      if (!reverse)
        ptr = 0;
    }
    if (ptr < 0 && reverse) {
      error("Out of boundary in looking for ibs0 in block.");
    }
  }

  if (start < 0 && !reverse) {start = (*posvec)[0];}
  if (start < 0 &&  reverse) {start = posvec->back();}

  bpos = (ptr > 0) ? ((ptr + 1) >> 3) : 0;
  k = bpos;

  // Process the block closest to / containing the starting pos
  uint8_t byte = ( iRR[k] ^ jRR[k] ) & ( iAA[k] ^ jAA[k] );
  if (byte) {
    if (reverse) {
      for (int32_t bit=0; bit<8 && k*8+bit<bmatRR->ncol; ++bit) {
        int32_t pt = k*8+(7-bit);
        if ( ((byte >> bit) & 0x01) && ((*posvec)[pt] <= start) ) {
          return (*posvec)[pt];
        }
      }
    }
    else {
      for (int32_t bit=0; bit<8 && k*8+bit<bmatRR->ncol; ++bit) {
        int32_t pt = k*8+bit;
        if ( ((byte >> (7-bit)) & 0x01) && ((*posvec)[pt] >= start) ) {
          return (*posvec)[pt];
        }
      }
    }
  }

  if (reverse && bpos == 0) {
    return -1;
  }
  if (!reverse && bpos >= bmatRR->nbytes_col - 1) {
    return -1;
  }
  if (reverse && bpos > 0) {
    k = bpos - 1;
    while(k >= 0) {
      if ( ( iRR[k] ^ jRR[k] ) & ( iAA[k] ^ jAA[k] ) ) { // IBS0 exists
        break;
      }
      k--;
    }
    if (k == -1)
      return -1;
  }
  if ((!reverse) && bpos < bmatRR->ncol - 1) {
    k=bpos + 1;
    while(k < bmatRR->nbytes_col) {
      if ( ( iRR[k] ^ jRR[k] ) & ( iAA[k] ^ jAA[k] ) ) { // IBS0 exists
        break;
      }
      k++;
    }
    if (k == bmatRR->nbytes_col)
      return -1;
  }

  // Find the position of IBS0
  byte = ( iRR[k] ^ jRR[k] ) & ( iAA[k] ^ jAA[k] );
  if (byte) {
    if (reverse) {
      for (int32_t bit=0; bit<8 && k*8+bit<bmatRR->ncol; ++bit) {
        int32_t pt = k*8+(7-bit);
        if ( ((byte >> bit) & 0x01) && ((*posvec)[pt] <= start) ) {
          return (*posvec)[pt];
        }
      }
    }
    else {
      for (int32_t bit=0; bit<8 && k*8+bit<bmatRR->ncol; ++bit) {
        int32_t pt = k*8+bit;
        if ( ((byte >> (7-bit)) & 0x01) && ((*posvec)[pt] >= start) ) {
          return (*posvec)[pt];
        }
      }
    }
  }
  warning("Something went wrong in IBS0inOneBlock. %d, %d, %d, %d (%d), %u", start, reverse, bpos, k, bmatRR->nbytes_col, byte);
  return -2;
}

// Read a block of genotypes
// for each variant with at least min_hom_gts homozygotes of both alleles
// store its position and homozygoes in compact matrices
int32_t ReadIBS0Block(std::string& reg, std::string& inVcf,
                   std::vector<bitmatrix*>& bmatRR_que,
                   std::vector<bitmatrix*>& bmatAA_que,
                   std::vector<std::vector<int32_t>* >& posvec_que,
                   int32_t min_hom_gts, bool add_to_front) {

  std::vector<GenomeInterval> ibs0interval;
  parse_intervals(ibs0interval, "", reg);
  BCFOrderedReader ibs0odr(inVcf, ibs0interval);
  int32_t nsamples = bcf_hdr_nsamples(ibs0odr.hdr);
  bcf1_t* ibs0iv = bcf_init();
  bitmatrix *bmatRR = new bitmatrix(nsamples, 1024);
  bitmatrix *bmatAA = new bitmatrix(nsamples, 1024);
  std::vector<int32_t> *posvec = new std::vector<int32_t>;
  uint8_t* gtRR = (uint8_t*)calloc(nsamples, sizeof(uint8_t));
  uint8_t* gtAA = (uint8_t*)calloc(nsamples, sizeof(uint8_t));
  int32_t* p_gt = NULL;
  int32_t  n_gt = 0;
  for(int32_t k=0; ibs0odr.read(ibs0iv); ++k) {  // read markers
    if ( bcf_get_genotypes(ibs0odr.hdr, ibs0iv, &p_gt, &n_gt) < 0 ) {
      notice("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(ibs0odr.hdr, ibs0iv->rid), ibs0iv->pos+1);
      return 0;
    }
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
    return 0;
  }

  if (add_to_front) {
    bmatRR_que.insert(bmatRR_que.begin(), bmatRR);
    bmatAA_que.insert(bmatAA_que.begin(), bmatAA);
    posvec_que.insert(posvec_que.begin(), posvec);
  } else {
    bmatRR_que.push_back(bmatRR);
    bmatAA_que.push_back(bmatAA);
    posvec_que.push_back(posvec);
  }
  return bmatRR->ncol;

}



IBS0lookup::IBS0lookup(const std::string &_inVcf, const std::string &_reg,
                       const bp2cmMap &_pgmap, double _margin,
                       int32_t _ck, int32_t _mh) :
  inVcf(_inVcf), reg(_reg), pgmap(_pgmap), margin_cM(_margin),
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
  HalfArmBound(start, end, leftend, rightend);

  ibs_st = mid - chunksize/2;
  ibs_ed = ibs_st + chunksize - 1;
  // Find the starting point for adding ibs0 lookup blocks
  while ( ibs_st > leftend && pgmap.bp2cm(start) - pgmap.bp2cm(ibs_st) < margin_cM ) {
    ibs_ed = ibs_st - 1;
    ibs_st -= chunksize;
  }
  if (ibs_st <= leftend) {
    reached_leftend = true;
    ibs_st = leftend;
  }

  // Moving forward, read common variants in chunks
  while (ibs_ed < rightend && pgmap.bp2cm(ibs_st)-pgmap.bp2cm(end) < margin_cM) {
    reg = chrom + ":" + std::to_string(ibs_st) + "-" + std::to_string(ibs_ed);
    if (ReadIBS0Block(reg, inVcf, bmatRR_que, bmatAA_que, posvec_que, min_hom_gts, 0) > 0) {
      start_que.push_back((*posvec_que.back())[0]);
    }
    ibs_st = ibs_ed + 1;
    ibs_ed +=chunksize;
  } // Finish initialize ibs0 lookup blocks
  if (ibs_ed >= rightend) {
    reached_rightend = true;
  }
}

IBS0lookup::IBS0lookup(const std::string &_inVcf, const std::string &_reg,
                       const bp2cmMap &_pgmap, int32_t _margin,
                       int32_t _ck, int32_t _mh) :
  inVcf(_inVcf), reg(_reg), pgmap(_pgmap), margin_bp(_margin),
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
  notice("IBS0lookup initialized for with %s", reg.c_str());
  chrom = v[0];
  mid = (start + end) / 2;
  HalfArmBound(start, end, leftend, rightend);

  // To jump over centromere
  int32_t off_left = 0, off_right = 0;
  if (pgmap.overlap_centro(start,end) && pgmap.centromere_st > 0) {
    off_left = std::max(start - pgmap.centromere_st, 0);
    off_right = std::max(pgmap.centromere_ed - end, 0);
  }

  // Locate the left most point to start
  ibs_st = mid - chunksize/2;
  ibs_ed = ibs_st + chunksize;
  while ( ibs_st > leftend & start - ibs_st < margin_bp + off_left ) {
    ibs_ed = ibs_st - 1;
    ibs_st -= chunksize;
  }
  if (ibs_st <= leftend) {
    reached_leftend = true;
    ibs_st = leftend;
  }

  // Moving forward
  while (ibs_ed < rightend && ibs_st - end < margin_bp + off_right) {
    std::string breg = chrom + ":" + std::to_string(ibs_st) + "-" + std::to_string(ibs_ed);
    notice("IBS0lookup read chunk %s", breg.c_str());
    if (ReadIBS0Block(breg, inVcf, bmatRR_que, bmatAA_que, posvec_que, min_hom_gts, 0) > 0) {
      start_que.push_back((*posvec_que.back())[0]);
    }
    ibs_st = ibs_ed + 1;
    ibs_ed +=chunksize;
  } // Finish initialize ibs0 lookup blocks
  if (ibs_ed >= rightend) {
    reached_rightend = true;
  }
}

void IBS0lookup::Clear() {
  for (uint32_t it = 0; it < bmatAA_que.size(); ++it) {
    delete bmatAA_que[it];
    delete bmatRR_que[it];
    delete posvec_que[it];
  }
  bmatAA_que.clear();
  bmatRR_que.clear();
  posvec_que.clear();
  start_que.clear();
}

// Find IBS0 for a pair of individuals starting from pos
int32_t IBS0lookup::FindIBS0 (int32_t i, int32_t j, int32_t pos, bool reverse) {

  if ((reverse && (*posvec_que[0])[0] > pos) ||
      (!reverse && posvec_que.back()->back() < pos)) {
    return -1;
  }
  int32_t k = (int32_t) start_que.size() - 1;
  while (start_que[k] > pos && k > 0) {k--;} // Locate the block covering pos
  int32_t ibs0 = IBS0inOneBlock(bmatRR_que[k], bmatAA_que[k], posvec_que[k],
                                i, j, reverse, pos);
  if (reverse) {
    while (ibs0 < 0 && k > 0) {
      k--; // mover left
      ibs0 = IBS0inOneBlock(bmatRR_que[k], bmatAA_que[k], posvec_que[k],
                            i, j, reverse);
    }
    // if (ibs0 < 0 && start_que[0] )
    if (ibs0 > 0 && pgmap.centromere_ed < pos && pgmap.centromere_ed > ibs0) {
      ibs0 = pgmap.centromere_ed;
    }
  } else {
    while(ibs0 < 0 && k < (int32_t) start_que.size() - 1) {
      k++; // move right
      ibs0 = IBS0inOneBlock(bmatRR_que[k], bmatAA_que[k], posvec_que[k],
                            i, j, reverse);
    }
    if (ibs0 >= 0 && pgmap.centromere_st > pos && pgmap.centromere_st < ibs0) {
      ibs0 = pgmap.centromere_st;
    }
  }
  return ibs0;
}

// Update to delete & add blocks to cover a new region
int32_t IBS0lookup::Update(const std::string &_reg) {

  reg = _reg;
  int32_t start, end;
  int32_t ibs_st, ibs_ed;
  std::vector<std::string> v;
  std::string chrom;
  split(v, ":-", reg);
  if (!str2int32(v[1], start) || !str2int32(v[2], end)) {
    error("Invalid region for update ibs0 lookup blocks.");
  }
  HalfArmBound(start, end, leftend, rightend);
  chrom = v[0];

  if (by_bp) {
    int32_t margin = margin_bp;
    // Delete unnecessary blocks
    while (posvec_que.size() > 0 &&
           start - posvec_que[0]->back() > margin) {
      delete posvec_que[0]; posvec_que.erase(posvec_que.begin());
      delete bmatRR_que[0]; bmatRR_que.erase(bmatRR_que.begin());
      delete bmatAA_que[0]; bmatAA_que.erase(bmatAA_que.begin());
      start_que.erase(start_que.begin());
    }
    while (posvec_que.size() > 0 &&
           start_que.back() - end > margin) {
      delete posvec_que.back(); posvec_que.pop_back();
      delete bmatRR_que.back(); bmatRR_que.pop_back();
      delete bmatAA_que.back(); bmatAA_que.pop_back();
      start_que.pop_back();
    }

    // Add new ibs0 lookup blocks to tail
    if (posvec_que.size() == 0) {
      ibs_ed = start - 1;
    } else {
      ibs_ed = posvec_que.back()->back();
    }
    while (ibs_ed - end < margin && ibs_ed < rightend) {
      ibs_st = ibs_ed + 1;
      ibs_ed = ibs_st + chunksize;
      std::string breg = chrom + ":" + std::to_string(ibs_st) + "-" + std::to_string(ibs_ed);
      if (ReadIBS0Block(breg, inVcf, bmatRR_que, bmatAA_que, posvec_que, min_hom_gts) > 0) {
        start_que.push_back((*posvec_que.back())[0]);
      }
    }
    if (ibs_ed >= rightend) {
      reached_rightend = true;
    }
    // Add new ibs0 lookup blocks to head
    if (posvec_que.size() == 0) {
      ibs_st = end + 1;
    } else {
      ibs_st = start_que[0];
    }
    while (start - ibs_st < margin && ibs_st > leftend) {
      ibs_ed = ibs_st - 1;
      ibs_st = ibs_ed - chunksize;
      std::string breg = chrom + ":" + std::to_string(ibs_st) + "-" + std::to_string(ibs_ed);
      if (ReadIBS0Block(breg, inVcf, bmatRR_que, bmatAA_que, posvec_que, min_hom_gts,1) > 0) {
        start_que.push_back((*posvec_que[0])[0]);
      }
    }
    if (ibs_st <= leftend) {
      reached_leftend = true;
    }
  }

  else {
    double margin = margin_cM;
    // Delete unnecessary blocks
    while (posvec_que.size() > 0 &&
           pgmap.bp2cm(start) - pgmap.bp2cm(posvec_que[0]->back()) > margin) {
      delete posvec_que[0]; posvec_que.erase(posvec_que.begin());
      delete bmatRR_que[0]; bmatRR_que.erase(bmatRR_que.begin());
      delete bmatAA_que[0]; bmatAA_que.erase(bmatAA_que.begin());
      start_que.erase(start_que.begin());
    }
    while (posvec_que.size() > 0 &&
           pgmap.bp2cm((*posvec_que.back())[0]) - pgmap.bp2cm(end) > margin) {
      delete posvec_que.back(); posvec_que.pop_back();
      delete bmatRR_que.back(); bmatRR_que.pop_back();
      delete bmatAA_que.back(); bmatAA_que.pop_back();
      start_que.pop_back();
    }

    // Add new ibs0 lookup blocks to tail
    if (posvec_que.size() == 0) {
      ibs_ed = start - 1;
    } else {
      ibs_ed = posvec_que.back()->back();
    }
    while (ibs_ed < rightend && pgmap.bp2cm(ibs_ed) - pgmap.bp2cm(end) < margin) {
      ibs_st = ibs_ed + 1;
      ibs_ed = ibs_st + chunksize;
      std::string breg = chrom + ":" + std::to_string(ibs_st) + "-" + std::to_string(ibs_ed);
      if (ReadIBS0Block(breg, inVcf, bmatRR_que, bmatAA_que, posvec_que, min_hom_gts) > 0){
        start_que.push_back((*posvec_que.back())[0]);
      }
    }
    if (ibs_ed >= rightend) {
      reached_rightend = true;
    }
    // Add new ibs0 lookup blocks to head
    if (posvec_que.size() == 0) {
      ibs_st = end + 1;
    } else {
      ibs_st = start_que[0];
    }
    while (ibs_st > leftend && pgmap.bp2cm(start) - pgmap.bp2cm(ibs_st) < margin) {
      ibs_ed = ibs_st - 1;
      ibs_st = ibs_ed - chunksize;
      std::string breg = chrom + ":" + std::to_string(ibs_st) + "-" + std::to_string(ibs_ed);
      if (ReadIBS0Block(breg, inVcf, bmatRR_que, bmatAA_que, posvec_que, min_hom_gts,1) > 0){
        start_que.push_back((*posvec_que[0])[0]);
      }
    }
    if (ibs_st <= leftend) {
      reached_leftend = true;
    }
  }
  return (int32_t) start_que.size();
}


// Update to include blocks for exactly the specified region
int32_t IBS0lookup::Update_Fixed(std::string &_reg) {
  int32_t start, end;
  std::vector<std::string> v;
  split(v, ":-", _reg);
  if (!str2int32(v[1], start) || !str2int32(v[2], end)) {
    error("Invalid region for initialize ibs0 lookup blocks. %s", _reg.c_str());
  }
  if (start >= end) {
    error("Invalid region for initialize ibs0 lookup blocks. %s", _reg.c_str());
  }
  Clear();
  if (start <= leftend) {
    reached_leftend = true;
  }
  if (end >= rightend) {
    reached_rightend = true;
  }
  if (start > rightend || end < leftend) {
    return -1;
  }
  if (ReadIBS0Block(_reg, inVcf, bmatRR_que, bmatAA_que, posvec_que, min_hom_gts) > 0) {
    start_que.push_back((*posvec_que.back())[0]);
  } else {
    notice("IBS0lookup::Update_Fixed Empty region %s", _reg.c_str());
  }
  return (int32_t) start_que.size();
}

void IBS0lookup::HalfArmBound(int32_t start, int32_t end, int32_t &leftend, int32_t &rightend) {
  if (start < pgmap.centromere_ed) {
    leftend = 0;
  } else {
    leftend = pgmap.centromere_ed;
  }
  if (end > pgmap.centromere_st) {
    rightend = pgmap.maxpos;
  } else {
    rightend = pgmap.centromere_st;
  }
}
