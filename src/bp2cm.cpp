#include "bp2cm.h"

bp2cmMap::bp2cmMap(const std::string &inMap, const char* sep,
                   int32_t cst, int32_t ced,
                   int32_t bp_col, int32_t cm_col,
                   int32_t _binsize) : binsize(_binsize) {

  std::map<int32_t, double> lookup;
  std::map<int32_t, std::vector<int32_t> > poshash;

  std::ifstream mfile (inMap);
  std::string line;
  std::vector<std::string> words;
  int32_t prepos = 0;
  double precm = 0.0;

  if (mfile.is_open()) {
    while(std::getline(mfile, line)) {
      split(words, sep, line);
      lookup[std::stoi(words[bp_col])] = std::stod(words[cm_col]);
    }
  }
  minpos = lookup.begin()->first;
  maxpos = std::stoi(words[bp_col]);
  maxcm  = std::stod(words[cm_col]);
// std::cout << minpos <<'\t' << maxpos << '\t' << maxcm << '\n';
  int32_t kb = 0;
  while (kb <= minpos / binsize) {
    poshash[kb] = std::vector<int32_t>{minpos};
    kb++;
  }

  for (auto it = lookup.begin(); it != lookup.end(); ++it) {
    kb = it->first / binsize;
    if ( poshash.find(kb) != poshash.end() ) {
      (poshash[kb]).push_back(it->first);
    } else {
      poshash[kb] =  std::vector<int32_t>{it->first};
    }
    if (cst < 0) {
      if (prepos > minpos && it->first - prepos > 200000 && it->second-precm < 1e-8) { // chr1-21
        cst = prepos; ced = it->first;
      } else if (minpos > 1e06) { // chr22
        cst = 10864561;
        ced = 12915808;
      } else {
        prepos = it->first;
        precm  = it->second;
      }
    }
  }
std::cout << "Centromere: " << cst <<'\t' <<ced << '\n';

  int32_t MAX = maxpos + binsize;
  for (kb = 0; kb < MAX/binsize; kb++) {
    if (poshash.find(kb) == poshash.end()) {
      if (kb*binsize > cst && kb*binsize < ced) {
        poshash[kb] = (kb*binsize-cst > ced-kb*binsize) ? std::vector<int32_t>{ced}:std::vector<int32_t>{cst};
      }
      else {
        std::vector<int32_t> rec;
        int32_t i = kb-1;
        while (poshash.find(i) == poshash.end() && i > 0) {
          i--;
        }
        int32_t j = kb+1;
        while (poshash.find(j) == poshash.end() && j < MAX/binsize) {
          j++;
        }
        // poshash[kb][0] = (kb-i < j-kb) ? poshash[i].back():poshash[j][0];
        if (kb-i == j-kb) {
          rec.push_back(poshash[i].back());
          rec.push_back(poshash[j][0]);
        } else if (kb-i < j-kb) {
          rec.push_back(poshash[i].back());
        } else {
          rec.push_back(poshash[j][0]);
        }
        poshash[kb] = rec;
      }
    }
  }
  for (kb = 0; kb < MAX/binsize; kb++) {
    if ( (kb*binsize < cst || kb*binsize > ced) && kb > 0) {
      if (kb*binsize - poshash[kb-1].back() > 0 && kb*binsize - poshash[kb-1].back() < poshash[kb][0] - kb*binsize ) {
        poshash[kb].push_back(poshash[kb-1].back());
      }
      if (kb*binsize - poshash[kb-1].back() > poshash[kb][0] - kb*binsize) {
        poshash[kb-1].push_back(poshash[kb][0]);
      }
    }
  }

  // Transport poshash to include genetic distance.
  for (auto & v : poshash) {
    std::vector<std::pair<int32_t, double> > rec;
    bphash[v.first] = rec;
    std::sort(v.second.begin(), v.second.end());
    for (auto & w : v.second) {
      bphash[v.first].push_back(std::pair<int32_t, double>(w, lookup[w]));
    }
  }

  centromere_st = cst;
  centromere_ed = ced;
}

double bp2cmMap::bp2cm(int32_t pos) {
  int32_t kb = pos/binsize;
  auto & v = bphash[kb];   // A vector of pairs
  int32_t mindist = (int32_t) 1e9;
  if (pos < minpos) {
    return 0.0;
  }
  if (pos > maxpos) {
    return maxcm;
  }
  double cm = 0.0;
  for (auto & rec : v) {
    if (std::abs(rec.first-pos) < mindist) {
      cm = rec.second;
      mindist = std::abs(rec.first-pos);
    }
  }
  return cm;
}





