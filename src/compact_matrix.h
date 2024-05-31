#ifndef __COMPACT_MATRIX_H
#define __COMPACT_MATRIX_H

#include <stdint.h>
#include <stddef.h>
#include <cstdlib>
#include <cstring>
#include "byte_operation.h"

// row-wise stored 1-bit matrix
class bitmatrix {
public:
  int32_t nrow;
  int32_t ncol;
  int32_t nrow_alloc;
  int32_t nbytes_col; // internal values
  int32_t max_row_step;
  uint8_t* bytes;
  bitmatrix(int32_t _ncol, int32_t _nrow_alloc = 100, int32_t _max_row_step = 50000);
  ~bitmatrix() { if ( bytes != NULL ) free(bytes); }
  int32_t reserve(int32_t new_nrow_alloc = 0);
  int32_t add_row_ints(int32_t* intarray);
  int32_t add_row_bytes(uint8_t* bytearray);
  int32_t add_row_bits(uint8_t* bitarray);
  bool transpose(bool verbose = 0);
  inline uint8_t* get_row_bits(int32_t irow) { return bytes + irow*nbytes_col; }
  inline uint8_t get_byte_at(int32_t irow, int32_t icol) { return ( ( bytes[irow*nbytes_col + (icol >> 3)] >> (icol & 0x03) ) & 0x01 ); }
  //void head(int32_t r, int32_t c);
  void print(int32_t rbeg, int32_t rend, int32_t cbeg, int32_t cend);
  void print_int(int32_t rbeg, int32_t rend, int32_t cbeg, int32_t cend);

  int32_t inner_prod_and_bytes(int32_t i, int32_t j);

};

#endif // __COMPACT_MATRIX_H
