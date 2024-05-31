#ifndef __BYTE_PAIR_OPERATION_H
#define __BYTE_PAIR_OPERATION_H

#include <stdint.h>
#include <stddef.h>
#include <cstdlib>
#include <cstring>

class BytePairOperation {
public:
  uint8_t byte_pair_xor[65536];
  uint8_t byte_pair_and[65536];

  BytePairOperation() {
    int32_t i, j, k;
    for(i=0; i < 256; ++i) {
      for(j=0; j < 256; ++j) {
        uint8_t x = ((uint8_t)i ^ (uint8_t)j) & 0x0ff;
        uint8_t s = 0;
        for(k = 0; k < 8; ++k) { s += ((x >> k) & 0x01); }
        byte_pair_xor[i*256+j] = s;
        x = ((uint8_t)i & (uint8_t)j) & 0x0ff;
        s = 0;
        for(k = 0; k < 8; ++k) { s += ((x >> k) & 0x01); }
        byte_pair_and[i*256+j] = s;
      }
    }
  }

  inline int32_t sum_xor(uint8_t* x, uint8_t* y, int32_t nbytes) {
    int32_t sum = 0;
    for(int32_t i=0; i < nbytes; ++i) sum += byte_pair_xor[(x[i] << 8) + y[i]];
    return sum;
  }

  inline int32_t sum_and(uint8_t* x, uint8_t* y, int32_t nbytes) {
    int32_t sum = 0;
    for(int32_t i=0; i < nbytes; ++i) sum += byte_pair_and[(x[i] << 8) + y[i]];
    return sum;
  }

  inline bool positive_xor(uint8_t* x, uint8_t* y, int32_t nbytes) {
    for(int32_t i=0; i < nbytes; ++i) {
      if ( byte_pair_xor[(x[i] << 8) + y[i]] ) return true;
    }
    return false;
  }

  inline bool positive_and(uint8_t* x, uint8_t* y, int32_t nbytes) {
    for(int32_t i=0; i < nbytes; ++i) {
      if ( byte_pair_and[(x[i] << 8) + y[i]] ) return true;
    }
    return false;
  }

};

extern BytePairOperation byte_pair_op;

#endif
