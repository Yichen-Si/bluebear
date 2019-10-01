#ifndef BPCM_H
#define BPCM_H

#include "utils.h"

#include <fstream>
#include <algorithm>

class bp2cmMap
{
  public:
    int32_t binsize;
    int32_t centromere_st, centromere_ed;
    int32_t minpos, maxpos;
    std::map<int32_t, std::vector<std::pair<int32_t, double> > > bphash;
    double maxcm;


  bp2cmMap(const std::string &inMap, const char* sep,
           int32_t cst = -1, int32_t ced = -1,
           int32_t bp_col = 3, int32_t cm_col = 2,
           int32_t _binsize=1000);

  double bp2cm(int32_t pos);

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
