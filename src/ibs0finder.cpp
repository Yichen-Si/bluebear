#include "ibs0.h"
#include "ibs0finder.h"

IBS0Finder::IBS0Finder(const std::string &_inVcf, const std::string &_reg,
                       gMap &_pgmap, double _margin,
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
    if (ReadIBS0Block(reg, inVcf, bmatRR_que, bmatAA_que, posvec_que, min_hom_gts, 0) > 0) {
      start_que.push_back((*posvec_que.back())[0]);
    }
    ibs_st = ibs_ed + 1;
    ibs_ed +=chunksize;
  } // Finish initialize ibs0 lookup blocks
}

IBS0Finder::IBS0Finder(const std::string &_inVcf, const std::string &_reg,
                       gMap &_pgmap, int32_t _margin,
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
    if (ReadIBS0Block(reg, inVcf, bmatRR_que, bmatAA_que, posvec_que, min_hom_gts, 0) > 0) {
      start_que.push_back((*posvec_que.back())[0]);
    }
    ibs_st = ibs_ed + 1;
    ibs_ed +=chunksize;
  } // Finish initialize ibs0 lookup blocks
}

// Find IBS0 for a pair of individuals starting from pos
int32_t IBS0Finder::FindIBS0 (int32_t i, int32_t j, int32_t pos, bool reverse) {

  if ((reverse && (*posvec_que[0])[0] > pos) ||
      (!reverse && posvec_que.back()->back() < pos)) {
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
int32_t IBS0Finder::Update(const std::string &_reg) {

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
      start_que.erase(start_que.begin());
    }
    while (start_que.back() - end > margin &&
           posvec_que.size() > 0) {
      delete posvec_que.back(); posvec_que.pop_back();
      delete bmatRR_que.back(); bmatRR_que.pop_back();
      delete bmatAA_que.back(); bmatAA_que.pop_back();
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
      if (ReadIBS0Block(reg, inVcf, bmatRR_que, bmatAA_que, posvec_que, min_hom_gts) > 0){
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
      if (ReadIBS0Block(reg, inVcf, bmatRR_que, bmatAA_que, posvec_que, min_hom_gts,1) > 0){
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
      start_que.erase(start_que.begin());
    }
    while (pgmap.bp2cm((*posvec_que.back())[0]) - pgmap.bp2cm(end) > margin &&
           posvec_que.size() > 0) {
      delete posvec_que.back(); posvec_que.pop_back();
      delete bmatRR_que.back(); bmatRR_que.pop_back();
      delete bmatAA_que.back(); bmatAA_que.pop_back();
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
      if (ReadIBS0Block(reg, inVcf, bmatRR_que, bmatAA_que, posvec_que, min_hom_gts) > 0){
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
      if (ReadIBS0Block(reg, inVcf, bmatRR_que, bmatAA_que, posvec_que, min_hom_gts,1) > 0){
        start_que.push_back((*posvec_que[0])[0]);
      }
    }
  }

  return (int32_t) start_que.size();

}


// Update to include blocks for exactly the region
int32_t IBS0Finder::Update_Fixed(const std::string &_reg) {
  reg = _reg;
  Clear();
  if (ReadIBS0Block(reg, inVcf, bmatRR_que, bmatAA_que, posvec_que, min_hom_gts) > 0) {
    start_que.push_back((*posvec_que.back())[0]);
  } else {
    notice("IBS0Finder::Update_Fixed Empty region");
  }
  return (int32_t) start_que.size();
}
