#include "ibs0.h"

int32_t IBS0inOneBlock(bitmatrix* bmatRR, bitmatrix* bmatAA,
                       int32_t i, int32_t j,
                       bool reverse, int32_t start) {
  uint8_t* iRR = bmatRR->get_row_bits(i);
  uint8_t* iAA = bmatAA->get_row_bits(i);
  uint8_t* jRR = bmatRR->get_row_bits(j);
  uint8_t* jAA = bmatAA->get_row_bits(j);
  int32_t k;
  if (reverse) {
    start = (start < 0) ? (bmatRR->nbytes_col-1) : start;
    for(k=start; k >= 0; --k) {
      if ( ( iRR[k] ^ jRR[k] ) & ( iAA[k] ^ jAA[k] ) ) { // IBS0 exists
        break;
      }
    }
    return std::min(k * 8 + 7, bmatRR->ncol-1);
  } else {
    start = (start < 0) ? 0 : start;
    for(k=start; k < bmatRR->nbytes_col; ++k) {
      if ( ( iRR[k] ^ jRR[k] ) & ( iAA[k] ^ jAA[k] ) ) { // IBS0 exists
        break;
      }
    }
    if (k == bmatRR->nbytes_col)
      return -1;
    return k * 8;
  }
}

int32_t ReadIBS0Block(std::string& reg, std::string& inVcf,
                   std::vector<bitmatrix*>& bmatRR_que,
                   std::vector<bitmatrix*>& bmatAA_que,
                   std::vector<std::vector<int32_t>* >& posvec_que,
                   int32_t min_hom_gts, bool add_to_front) {

    std::vector<GenomeInterval> ibs0interval;
    parse_intervals(ibs0interval, "", reg);
    BCFOrderedReader ibs0odr(inVcf, ibs0interval);
    int32_t nsamples = bcf_hdr_nsamples(ibs0odr.hdr);
    bcf1_t* ibs0iv = bcf_init();
    bitmatrix *bmatRR = new bitmatrix(nsamples);
    bitmatrix *bmatAA = new bitmatrix(nsamples);
    std::vector<int32_t> *posvec = new std::vector<int32_t>;
    uint8_t* gtRR = (uint8_t*)calloc(nsamples, sizeof(uint8_t));
    uint8_t* gtAA = (uint8_t*)calloc(nsamples, sizeof(uint8_t));
    int32_t* p_gt = NULL;
    int32_t  n_gt = 0;
    for(int32_t k=0; ibs0odr.read(ibs0iv); ++k) {  // read markers
      if ( bcf_get_genotypes(ibs0odr.hdr, ibs0iv, &p_gt, &n_gt) < 0 ) {
        notice("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__FUNCTION__, bcf_hdr_id2name(ibs0odr.hdr, ibs0iv->rid), ibs0iv->pos+1);
        return 0;
      }
      memset(gtRR, 0, nsamples);
      memset(gtAA, 0, nsamples);
      int32_t gcs[3] = {0,0,0};
      for(int32_t i=0; i < nsamples; ++i) {
        int32_t g1 = p_gt[2*i];
        int32_t g2 = p_gt[2*i+1];
        int32_t geno;
        if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
          //geno = 0;
        } else {
          geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0);
          if ( geno == 0 )    { gtRR[i] = 1; }
          else if ( geno == 2 ) { gtAA[i] = 1; }
          ++gcs[geno];
        }
      }
      if (( gcs[0] >= min_hom_gts ) && ( gcs[2] >= min_hom_gts )) {
        bmatRR->add_row_bytes(gtRR);
        bmatAA->add_row_bytes(gtAA);
        posvec->push_back(ibs0iv->pos+1);
      }
    }
    free(gtRR);
    free(gtAA);
    bmatRR->transpose();
    bmatAA->transpose();
    if (bmatRR->ncol == 0) {
      return 0;
    }

    if (add_to_front) {
      bmatRR_que.insert(bmatRR_que.begin(), bmatRR);
      bmatAA_que.insert(bmatAA_que.begin(), bmatAA);
      posvec_que.insert(posvec_que.begin(), posvec);
    } else {
      bmatRR_que.push_back(bmatRR);
      bmatAA_que.push_back(bmatAA);
      posvec_que.push_back(posvec);
    }
    return bmatRR->ncol;

  }
