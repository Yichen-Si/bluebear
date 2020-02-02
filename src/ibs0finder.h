
#ifndef IBS0Finder_H
#define IBS0Finder_H

#include "cramore.h"
#include "bcf_ordered_reader.h"
#include "compact_matrix.h"
#include "utils.h"
#include "gmap.h"

class IBS0Finder
{
public:
  std::string inVcf, reg;
  std::vector<bitmatrix*> bmatRR_que, bmatAA_que;
  std::vector<std::vector<int32_t>* > posvec_que;
  std::vector<int32_t> start_que;
  gMap pgmap;
  double margin_cM;
  int32_t margin_bp;
  int32_t chunksize, min_hom_gts;
  bool by_bp = 0;

  IBS0Finder(const std::string &_inVcf, const std::string &_reg,
             gMap &_pgmap, double _margin,
             int32_t _ck = 1000000, int32_t _mh = 1);
  IBS0Finder(const std::string &_inVcf, const std::string &_reg,
             gMap &_pgmap, int32_t _margin,
             int32_t _ck = 1000000, int32_t _mh = 1);

  // Find IBS0 for a pair of individuals to one direction starting from pos
  int32_t FindIBS0 (int32_t i, int32_t j, int32_t pos, bool reverse = 0);

  // Update to delete & add blocks to cover a new region
  int32_t Update (const std::string &_reg);
  int32_t Update_Fixed (const std::string &_reg);
  void Clear() {
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

  ~IBS0Finder() {Clear();}

};


#endif
