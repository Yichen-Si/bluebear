#ifndef BPCM_H
#define BPCM_H

#include "utils.h"

#include <fstream>
#include <algorithm>

struct mapEntry
{
  int32_t pos;
  double rate, cm;
  mapEntry(int32_t _p, double _r, double _c) : pos(_p), rate(_r), cm(_c) {}
};

class gMap
{
  public:
    int32_t binsize;
    int32_t centromere_st, centromere_ed;
    int32_t minpos, maxpos;
    std::map<int32_t, std::vector<mapEntry*> > bphash;
    double maxcm, ctrcm;

  gMap(const std::string &inMap, const char* sep,
       std::string chrom, int32_t cst = -1, int32_t ced = -1,
       int32_t bp_col = 1, int32_t cm_col = 3, int32_t rate_col = 2,
       int32_t _binsize=1000);

  double bp2cm(int32_t pos);
  double bpinterval2cm(int32_t st, int32_t ed);

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
