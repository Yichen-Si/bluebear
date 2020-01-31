#include "gmap.h"

gMap::gMap(const std::string &inMap, const char* sep,
           std::string chrom, int32_t cst, int32_t ced,
           int32_t bp_col, int32_t cm_col, int32_t rate_col,
           int32_t _binsize) : binsize(_binsize) {

  std::ifstream mfile (inMap);
  std::string line;
  std::vector<std::string> words;
  int32_t prepos = 0;
  double precm = 0.0;

  if (mfile.is_open()) {
    while(std::getline(mfile, line)) {
      split(words, sep, line);
      if (words[0] != chrom) {continue;}
      prepos = std::stoi(words[bp_col]);
      mapEntry *rec = new mapEntry(prepos, std::stod(words[rate_col]), std::stod(words[cm_col]));
      bphash[prepos/binsize].push_back(rec);
    }
  }
  minpos = bphash.begin()->second[0]->pos;
  maxpos = bphash.rbegin()->second.back()->pos;
  maxcm  = bphash.rbegin()->second.back()->cm;
  prepos = 0;
  if (cst == -1 && ced == -1) { // Find centromere
    for (auto it = bphash.begin(); it != bphash.end(); ++it) {
      for (auto pt : it->second) {
        if (cst < 0) {
          if (prepos > minpos && pt->pos - prepos > 200000 && pt->cm - precm < 1e-8) { // chr1-21
            cst = prepos; ced = pt->pos; ctrcm = pt->cm;
          } else if (minpos > 1e06) { // chr22
            cst = 10864561;
            ced = 15344143;
            ctrcm = 0.01;
          } else {
            prepos = pt->pos;
            precm  = pt->cm;
          }
        }
      }
    }
  }
  centromere_st = cst;
  centromere_ed = ced;
}


double gMap::bp2cm(int32_t pos) {

  if (pos <= minpos) {return 0.0;}
  if (pos >= maxpos) {return maxcm;}
  if (in_centro(pos)){return ctrcm;}

  // Find key
  int32_t key = pos/binsize;
  auto rec = bphash.find(key);

  // Find the closest position upstream
  int32_t st_bp = 0;
  double  st_cm = 0.0;
  if (rec != bphash.end()) {
    for (auto it = rec->second.rbegin(); it != rec->second.rend(); ++it) {
      if ((*it)->pos <= pos) {
        st_bp = (*it)->pos;
        st_cm = (*it)->cm + ((*it)->rate)*((pos-st_bp)*1e-6);
        break;
      }
    }
  }
  while(rec == bphash.end() && st_bp == 0) {
    key--;
    rec = bphash.find(key);
  }
  st_bp = rec->second.back()->pos;
  st_cm = rec->second.back()->cm;
  st_cm += (rec->second.back()->rate)*((pos-st_bp)*1e-6);

  return st_cm;
}


double gMap::bpinterval2cm(int32_t st, int32_t ed) {

  int32_t st_bp = -1, ed_bp = -1;
  double  st_cm = -1.0, ed_cm = -1.0;
  int32_t key;

  if (in_centro(st, ed)) {return 0.0;}
  if (ed <= minpos) {return 0.0;}
  if (st >= maxpos) {return 0.0;}
  if (st <= minpos) {st_bp = minpos; st_cm = 0.0;}
  if (ed >= maxpos) {ed_bp = maxpos; ed_cm = maxcm;}
  if (in_centro(st)) {st_cm = ctrcm;}
  if (in_centro(ed)) {ed_cm = ctrcm;}

  key = st/binsize;
  auto rec = bphash.find(key);
  if (st_cm < 0.0) {
    // Find the genetic position of st
    if (rec != bphash.end()) {
      for (auto it = rec->second.rbegin(); it != rec->second.rend(); ++it) {
        if ((*it)->pos <= st) {
          st_bp = (*it)->pos;
          st_cm = (*it)->cm + ((*it)->rate)*((st-st_bp)*1e-6);
          break;
        }
      }
    }
    while(rec == bphash.end() && st_bp == 0) {
      key--;
      rec = bphash.find(key);
    }
    st_bp = rec->second.back()->pos;
    st_cm = rec->second.back()->cm;
    st_cm -= (rec->second.back()->rate)*((st-st_bp)*1e-6);
  }

  if (ed_cm < 0.0) {
    // Find the genetic position of ed
    key = ed/binsize;
    rec = bphash.find(key);
    if (rec != bphash.end()) {
      for (auto it = rec->second.begin(); it != rec->second.end(); ++it) {
        if ((*it)->pos >= ed) {
          ed_bp = (*it)->pos;
          ed_cm = (*it)->cm - ((*it)->rate)*((ed_bp-ed)*1e-6);
          break;
        }
      }
    } else {
      while(rec == bphash.end() && ed_cm < 0.0) {
        key++;
        rec = bphash.find(key);
      }
      ed_bp = rec->second[0]->pos;
      ed_cm = rec->second[0]->cm;
      ed_cm -= (rec->second[0]->rate)*((ed_bp-ed)*1e-6);
    }
  }

  return ed_cm - st_cm;
}
















