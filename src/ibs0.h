
#ifndef IBS0_H
#define IBS0_H

#include "cramore.h"
#include "bcf_ordered_reader.h"
#include "compact_matrix.h"
#include "bp2cm.h"
#include "utils.h"

int32_t IBS0inOneBlock(bitmatrix* bmatRR, bitmatrix* bmatAA,
                       std::vector<int32_t>* posvec,
                       int32_t i, int32_t j, bool reverse = 0,
                       int32_t start = -1, int32_t _ptr = -1);

int32_t ReadIBS0Block(std::string& reg, std::string& inVcf,
                   std::vector<bitmatrix*>& bmatRR_que,
                   std::vector<bitmatrix*>& bmatAA_que,
                   std::vector<std::vector<int32_t>* >& posvec_que,
                   int32_t min_hom_gts=1, bool add_to_front=0);

class IBS0lookup
{
public:
  std::string inVcf, reg;
  std::vector<bitmatrix*> bmatRR_que, bmatAA_que;
  std::vector<std::vector<int32_t>* > posvec_que;
  std::vector<int32_t> start_que;
  bp2cmMap pgmap;
  double margin_cM;
  int32_t margin_bp;
  int32_t chunksize, min_hom_gts;
  bool by_bp = 0;

  IBS0lookup(const std::string &_inVcf, const std::string &_reg,
             bp2cmMap &_pgmap, double _margin,
             int32_t _ck = 1000000, int32_t _mh = 1);
  IBS0lookup(const std::string &_inVcf, const std::string &_reg,
             bp2cmMap &_pgmap, int32_t _margin,
             int32_t _ck = 1000000, int32_t _mh = 1);

  // Find IBS0 for a pair of individuals to both direction starting from pos
  int32_t FindIBS0 (int32_t i, int32_t j, int32_t pos, bool reverse = 0);

  // Update to delete & add blocks to cover a new region
  int32_t Update (const std::string &_reg);

  ~IBS0lookup() {
    for (uint32_t it = 0; it < start_que.size(); ++it) {
      delete bmatAA_que[it];
      delete bmatRR_que[it];
      delete posvec_que[it];
    }
  }

};

#endif
