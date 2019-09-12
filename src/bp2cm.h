#ifndef BPCM_H
#define BPCM_H

#include "utils.h"

#include <fstream>

class bp2cmMap
{
  public:
    std::map<int32_t, double> lookup;
    std::map<int32_t, std::vector<int32_t> > poshash;
    int32_t binsize;
    int32_t centromere_st, centromere_ed;
    int32_t minpos, maxpos;
    double maxcm;


  bp2cmMap(const std::string &inMap, const char* sep,
        int32_t bp_col = 3, int32_t cm_col = 2,
        int32_t _binsize=1000) : binsize(_binsize) {

    std::ifstream mfile (inMap);
    std::string line;
    std::vector<std::string> words;
    int32_t prepos = 0, cst = 0, ced = 0;

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
    while (kb <= lookup.begin()->first / binsize) {
      poshash[kb] = std::vector<int32_t>{lookup.begin()->first};
      kb++;
    }

    for (auto it = lookup.begin(); it != lookup.end(); ++it) {
      kb = it->first / binsize;
      if ( poshash.find(kb) != poshash.end() ) {
        (poshash[kb]).push_back(it->first);
      } else {
        poshash[kb] =  std::vector<int32_t>{it->first};
      }
      if (cst == 0) {
        if (it->first - prepos > 200000) {
          cst = prepos; ced = it->first;
        } else {
          prepos = it->first;
        }
      }
    }

    int32_t MAX = lookup.rbegin()->first;
    for (kb = 0; kb < MAX/binsize; kb++) {
      if (poshash.find(kb) == poshash.end()) {
        if (kb*binsize > cst && kb*binsize < ced) {
          poshash[kb] = (kb*binsize-cst > ced-kb*binsize) ? std::vector<int32_t>{ced}:std::vector<int32_t>{cst};
        }
        else {
          poshash[kb] = std::vector<int32_t>{0};
          int32_t i = kb-1;
          while (poshash.find(i) == poshash.end() && i > 0) {
            i--;
          }
          int32_t j = kb+1;
          while (poshash.find(j) == poshash.end() && j < MAX/binsize) {
            j++;
          }
          poshash[kb][0] = (kb-i < j-kb) ? poshash[i].back():poshash[j][0];
          if (kb-i == j-kb) {
            poshash[kb].push_back(poshash[i].back());
          }
        }
      }
      if ( (kb*binsize < cst || kb*binsize > ced) && kb > 0) {
        if (kb*binsize - poshash[kb-1].back() > 0 && kb*binsize - poshash[kb-1].back() < poshash[kb][0] - kb*binsize ) {
          poshash[kb].insert(poshash[kb].begin(), poshash[kb-1].back());
        }
        if (kb*binsize - poshash[kb-1].back() > poshash[kb][0] - kb*binsize) {
          poshash[kb-1].push_back(poshash[kb][0]);
        }
      }
    }
    centromere_st = cst;
    centromere_ed = ced;
  }

  double bp2cm(int32_t pos) {
    int32_t kb = pos/binsize;
    int32_t mindist = INT_MAX;
    if (pos < minpos) {
      return 0.0;
    }
    if (pos > maxpos) {
      return maxcm;
    }
    double cm = 0.0;
    for (auto&& rec : poshash[kb]) {
      if (std::abs(rec-pos) < mindist) {
        cm = lookup[rec];
        mindist = std::abs(rec-pos);
      }
    }
    return cm;
  }

  bool in_centro(int32_t pos) {
    return (pos > centromere_st && pos < centromere_ed);
  }

  bool in_centro(int32_t start, int32_t end) {
    return (start > centromere_st && end < centromere_ed);
  }

  int32_t overlap_centro(int32_t start, int32_t end) {
    return std::max(std::min(centromere_ed, end) - std::max(centromere_st, start), 0);
  }

};

#endif

