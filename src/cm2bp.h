#ifndef CMBP_H
#define CMBP_H

#include "utils.h"

#include <fstream>
#include <math.h>

class cm2bpMap
{
  public:
    int32_t scale;
    double binsize;
    std::map<int32_t, int32_t> lookup;
    int32_t minpos, maxpos;
    double maxcm;

  cm2bpMap(const std::string &inMap, const char* sep,
           int32_t bp_col = 3, int32_t cm_col = 2,
           double _binsize=0.1) : binsize(_binsize) {

    std::ifstream mfile (inMap);
    std::string line;
    std::vector<std::string> words;
    std::vector<std::vector<double> > lmap;

    // Read linkage map
    if (mfile.is_open()) {
      while(std::getline(mfile, line)) {
        split(words, sep, line);
        lmap.push_back(std::vector<double>{std::stod(words[cm_col]), std::stod(words[bp_col])} );
      }
    }
    minpos = (int32_t) lmap[0][1];
    maxpos = (int32_t) lmap.back()[1];
    maxcm  = lmap.back()[0];

    // Build lookup table
    scale = 1;
    double res = binsize * 0.1;
    int32_t pre = 0, k = 1;
    while (std::abs(binsize * scale - std::round(binsize * scale)) > res)
      scale *= 10;

    double tag = binsize * k;
    for (int32_t cur = 0; cur < (int32_t) lmap.size(); ++cur) {
      if (lmap[cur][0] < tag) {
        pre=cur;
        continue;
      }
      while (lmap[cur][0] > tag && k*binsize < lmap.back()[0]) {
        if (abs(lmap[cur][0] - tag) > abs(lmap[pre][0] - tag)) {
          lookup[(int32_t) k * binsize * scale] = (int32_t) lmap[pre][1];
        } else {
          lookup[(int32_t) k * binsize * scale] = (int32_t) lmap[cur][1];
        }
        k++;
        tag = binsize * k;
      }
      pre=cur;
    }
  }

  int32_t cm2bp(double cm) {
    if (cm > maxcm)
      return maxpos;
    return lookup[(int32_t) std::round(cm*scale)];
  }

};

#endif




