
#ifndef IBS0_H
#define IBS0_H

#include "cramore.h"
#include "bcf_ordered_reader.h"
#include "compact_matrix.h"

int32_t IBS0inOneBlock(bitmatrix* bmatRR, bitmatrix* bmatAA,
                       int32_t i, int32_t j,
                       bool reverse = 0, int32_t start = -1);

int32_t ReadIBS0Block(std::string& reg, std::string& inVcf,
                   std::vector<bitmatrix*>& bmatRR_que,
                   std::vector<bitmatrix*>& bmatAA_que,
                   std::vector<std::vector<int32_t>* >& posvec_que,
                   int32_t min_hom_gts=1, bool add_to_front=0);

#endif
