#include "bp2cm.h"

bp2cmMap::bp2cmMap(const std::string &inMap, const char* sep,
           std::string chrom, int32_t cst, int32_t ced,
           int32_t bp_col, int32_t cm_col, int32_t rate_col,
           int32_t _binsize) : binsize(_binsize) {

  std::ifstream mfile (inMap);
  std::string line;
  std::vector<std::string> words;
  int32_t prepos = 0, pos = 0, jump = 0, maxjump = 200000;
  double precm = 0.0, cm = 0.0, prerate = 0.0;

  if (!mfile.is_open()) {
    error("Genetic map file cannot be opened");
  }
  while(std::getline(mfile, line)) {
    split(words, sep, line);
    if (!chrom.empty() && words[0] != chrom) { // Need to filter for chromosome
      continue;
    }
    try {
      pos = std::stoi(words[bp_col]);
      cm  = std::stod(words[cm_col]);
      prerate = std::stod(words[rate_col]);
    }
    catch (const std::invalid_argument& ia) {
      continue;
    }
    jump = pos - prepos;
    if (pos > 1e6 && jump > maxjump) {
      cst = prepos; ced = pos; ctrcm = (precm+cm)/2.0;
      maxjump = jump;
    }
    precm  = cm;
    prepos = pos;
    mapEntry *rec = new mapEntry(pos, prerate, cm);
    bphash[prepos/binsize].push_back(rec);
  }

  // if (cst == 0) { // Temporary
  //   minpos = 500000;
  // } else {
  //   minpos = bphash.begin()->second[0]->pos;
  // }
  minpos = bphash.begin()->second[0]->pos;
  maxpos = bphash.rbegin()->second.back()->pos;
  maxcm  = bphash.rbegin()->second.back()->cm;
  centromere_st = cst;
  centromere_ed = ced;
}

void bp2cmMap::clear() {
  for (auto it = bphash.begin(); it != bphash.end(); ++it) {
    for (auto & pt : it->second) {
      delete pt; pt = NULL;
    }
  }
  notice("Cleared bp2cmMap object");
}

double bp2cmMap::bp2cm(int32_t pos) {

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


double bp2cmMap::bpinterval2cm(int32_t st, int32_t ed) {

  int32_t st_bp = -1, ed_bp = -1;
  double  st_cm = -1.0, ed_cm = -1.0;
  int32_t key;
  if(st > ed) {
    notice("bp2cmMap::bpinterval2cm: Invalid interval: %d,%d. End points are switched.",st,ed);
    int32_t tmp = ed;
    ed = st; st = tmp;
  }

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
    while(rec == bphash.end() && st_cm < 0.0) {
      key--;
      rec = bphash.find(key);
    }
    st_bp = rec->second.back()->pos;
    st_cm = rec->second.back()->cm;
    st_cm += (rec->second.back()->rate)*((st-st_bp)*1e-6);
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
    }
    while(rec == bphash.end() && ed_cm < 0.0) {
      key++;
      rec = bphash.find(key);
    }
    ed_bp = rec->second[0]->pos;
    ed_cm = rec->second[0]->cm;
    ed_cm -= (rec->second[0]->rate)*((ed_bp-ed)*1e-6);
  }

  if (ed_cm-st_cm < 0) { // Should never happen
    notice("bp2cmMap::bpinterval2cm: negtaive cm length %d,%.3f-%d,%.3f", st,st_cm,ed,ed_cm);
  }
  return ed_cm - st_cm;
}




// Replaced by new implement including recombination rate
// bp2cmMap::bp2cmMap(const std::string &inMap, const char* sep,
//                    int32_t cst, int32_t ced,
//                    int32_t bp_col, int32_t cm_col,
//                    int32_t _binsize) : binsize(_binsize) {

//   std::map<int32_t, double> lookup;
//   std::map<int32_t, std::vector<int32_t> > poshash;

//   std::ifstream mfile (inMap);
//   std::string line;
//   std::vector<std::string> words;
//   int32_t prepos = 0;
//   double precm = 0.0;

//   if (mfile.is_open()) {
//     while(std::getline(mfile, line)) {
//       split(words, sep, line);
//       lookup[std::stoi(words[bp_col])] = std::stod(words[cm_col]);
//     }
//   }
//   minpos = lookup.begin()->first;
//   maxpos = std::stoi(words[bp_col]);
//   maxcm  = std::stod(words[cm_col]);
// // std::cout << minpos <<'\t' << maxpos << '\t' << maxcm << '\n';
//   int32_t kb = 0;
//   while (kb <= minpos / binsize) {
//     poshash[kb] = std::vector<int32_t>{minpos};
//     kb++;
//   }

//   for (auto it = lookup.begin(); it != lookup.end(); ++it) {
//     kb = it->first / binsize;
//     if ( poshash.find(kb) != poshash.end() ) {
//       (poshash[kb]).push_back(it->first);
//     } else {
//       poshash[kb] =  std::vector<int32_t>{it->first};
//     }
//     if (cst < 0) {
//       if (prepos > minpos && it->first - prepos > 200000 && it->second-precm < 1e-8) { // chr1-21
//         cst = prepos; ced = it->first;
//       } else if (minpos > 1e06) { // chr22
//         cst = 10864561;
//         ced = 12915808;
//       } else {
//         prepos = it->first;
//         precm  = it->second;
//       }
//     }
//   }
// std::cout << "Centromere: " << cst <<'\t' <<ced << '\n';

//   int32_t MAX = maxpos + binsize;
//   for (kb = 0; kb < MAX/binsize; kb++) {
//     if (poshash.find(kb) == poshash.end()) {
//       if (kb*binsize > cst && kb*binsize < ced) {
//         poshash[kb] = (kb*binsize-cst > ced-kb*binsize) ? std::vector<int32_t>{ced}:std::vector<int32_t>{cst};
//       }
//       else {
//         std::vector<int32_t> rec;
//         int32_t i = kb-1;
//         while (poshash.find(i) == poshash.end() && i > 0) {
//           i--;
//         }
//         int32_t j = kb+1;
//         while (poshash.find(j) == poshash.end() && j < MAX/binsize) {
//           j++;
//         }
//         // poshash[kb][0] = (kb-i < j-kb) ? poshash[i].back():poshash[j][0];
//         if (kb-i == j-kb) {
//           rec.push_back(poshash[i].back());
//           rec.push_back(poshash[j][0]);
//         } else if (kb-i < j-kb) {
//           rec.push_back(poshash[i].back());
//         } else {
//           rec.push_back(poshash[j][0]);
//         }
//         poshash[kb] = rec;
//       }
//     }
//   }
//   for (kb = 0; kb < MAX/binsize; kb++) {
//     if ( (kb*binsize < cst || kb*binsize > ced) && kb > 0) {
//       if (kb*binsize - poshash[kb-1].back() > 0 && kb*binsize - poshash[kb-1].back() < poshash[kb][0] - kb*binsize ) {
//         poshash[kb].push_back(poshash[kb-1].back());
//       }
//       if (kb*binsize - poshash[kb-1].back() > poshash[kb][0] - kb*binsize) {
//         poshash[kb-1].push_back(poshash[kb][0]);
//       }
//     }
//   }

//   // Transport poshash to include genetic distance.
//   for (auto & v : poshash) {
//     std::vector<std::pair<int32_t, double> > rec;
//     bphash[v.first] = rec;
//     std::sort(v.second.begin(), v.second.end());
//     for (auto & w : v.second) {
//       bphash[v.first].push_back(std::pair<int32_t, double>(w, lookup[w]));
//     }
//   }

//   centromere_st = cst;
//   centromere_ed = ced;
// }

// double bp2cmMap::bp2cm(int32_t pos) {
//   int32_t kb = pos/binsize;
//   auto & v = bphash[kb];   // A vector of pairs
//   int32_t mindist = (int32_t) 1e9;
//   if (pos < minpos) {
//     return 0.0;
//   }
//   if (pos > maxpos) {
//     return maxcm;
//   }
//   double cm = 0.0;
//   for (auto & rec : v) {
//     if (std::abs(rec.first-pos) < mindist) {
//       cm = rec.second;
//       mindist = std::abs(rec.first-pos);
//     }
//   }
//   return cm;
// }

