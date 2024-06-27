#ifndef __BINARY_VARIANT_H
#define __BINARY_VARIANT_H

#include <set>
#include <tuple>
#include <iostream>
#include <vector>
#include <functional>
#include <unordered_map>
#include <random>
#include <algorithm>
#include <cstring>
#include <fstream>
#include "Error.h"
#include "byte_operation.h"

/**
 * Represent a binary variant in either sparse or dense format
 * Sparse: store indices of 1s
 * Dense: store a bit array
*/
struct BinaryVariant {
    bool isSparse;
    uint32_t nBits;  // Number of bits
    uint32_t nBytes; // Number of bytes
    uint32_t mOnes;  // Number of 1s in bitArray or indices.size()
    std::set<uint32_t> indices;  // always store the minory indices
    uint8_t* bitArray = nullptr; // fill from the most significant bit
    bool sign; //
    uint8_t mask;   // mask the last few bits in case nBits % 8 != 0

    BinaryVariant(bool _isSparse, uint32_t _n = 0, uint32_t _m = 0, bool _s = 0) : isSparse(_isSparse), nBits(_n), mOnes(_m), sign(_s) {
        nBytes = (nBits + 7) / 8;
        mask = 255 << (8 - nBits % 8);
    }
    BinaryVariant(std::set<uint32_t>& _indices, uint32_t _n = 0, bool _s = 0) : isSparse(true), nBits(_n), indices(_indices), sign(_s) {
        nBytes = (nBits + 7) / 8;
        mask = 255 << (8 - nBits % 8);
        mOnes = indices.size();
    }
    BinaryVariant(bool _isSparse, const int32_t* gt, uint32_t _n, bool zero = 0, bool _s = 0) : isSparse(_isSparse), nBits(_n), sign(_s) {
        nBytes = (nBits + 7) / 8;
        mask = 255 << (8 - nBits % 8);
        update_int<int32_t>(gt, zero);
    }
    BinaryVariant(bool _isSparse, const uint8_t* _bitArray, uint32_t _n, uint32_t _m = 0, bool zero = 0, bool _s = 0) : isSparse(_isSparse), nBits(_n), sign(_s) {
        nBytes = (nBits + 7) / 8;
        mask = 255 << (8 - nBits % 8);
        update_bit(_bitArray, _m, zero);
    }

    ~BinaryVariant() {
        if (bitArray != nullptr) {
            delete[] bitArray;
            bitArray = nullptr;
        }
    }

    // move constructor
    BinaryVariant(BinaryVariant&& other) noexcept : isSparse(other.isSparse), nBits(other.nBits), mOnes(other.mOnes), sign(other.sign) {
        nBytes = (nBits + 7) / 8;
        mask = 255 << (8 - nBits % 8);
        if (isSparse) {
            indices = std::move(other.indices);
        } else {
            bitArray = other.bitArray;
            other.bitArray = nullptr;
        }
    }

    // Assign operator
    BinaryVariant& operator=(BinaryVariant&& other) noexcept {
        if (this != &other) {
            if (bitArray != nullptr) {
                delete[] bitArray;
            }
            isSparse = other.isSparse;
            mOnes = other.mOnes;
            nBits = other.nBits;
            nBytes = (nBits + 7) / 8;
            sign = other.sign;
            if (isSparse) {
                indices = std::move(other.indices);
            } else {
                bitArray = other.bitArray;
                other.bitArray = nullptr;
            }
        }
        return *this;
    }

    // Currently only support comparison between the same format
    bool equal(const BinaryVariant& rhs) const {
        if (mOnes != rhs.mOnes || isSparse != rhs.isSparse) {
            return false;
        }
        if (isSparse && rhs.isSparse) {
            return indices == rhs.indices && sign == rhs.sign;
        }
        if (nBits != rhs.nBits) {
            return false;
        }
        if (sign != rhs.sign) { // prefer to never happen
            for (uint32_t i = 0; i < nBytes; ++i) {
                if (bitArray[i] != ~rhs.bitArray[i]) {
                    return false;
                }
            }
            return true;
        }
        return memcmp(bitArray, rhs.bitArray, nBytes) == 0;
    }

    void update_bit(const uint8_t* _bitArray, uint32_t _m = 0, bool _zero = 0) {
        if (bitArray != nullptr) {
            delete[] bitArray;
            bitArray = nullptr;
        }
        if (isSparse) {
            indices.clear();
            for (uint32_t i = 0; i < nBytes; ++i) {
                for (uint32_t j = 0; j < 8 && (i*8+j < nBits); ++j) {
                    if ((_bitArray[i] & (1 << (7 - j))) > 0 != _zero) {
                        indices.insert(i*8+j);
                    }
                }
            }
            mOnes = indices.size();
            return;
        }
        bitArray = new uint8_t[nBytes];
        if (_zero) {
            for (uint32_t i = 0; i < nBytes; ++i) {
                bitArray[i] = ~_bitArray[i];
            }
        } else {
            std::memcpy(bitArray, _bitArray, nBytes);
        }
        if (_m > 0) {
            mOnes = _m;
            return;
        }
        mOnes = 0; // if not specified, count 1s (after flipping if _zero)
        for (uint32_t i = 0; i < nBytes; ++i) {
            for (uint32_t j = 0; j < 8 && (i*8+j < nBits); ++j) {
                mOnes += (bitArray[i] & (1 << (7 - j))) != 0;
            }
        }
    }

    // Transform & store the input integer array in internal format
    // If sparse, store the indices !zero bits
    // If dense, store a bit array but flip the bits if !zero
    template <typename T>
    void update_int(const T* intArray, bool _zero = 0) {
        if (bitArray != nullptr) {
            delete[] bitArray;
            bitArray = nullptr;
        }
        indices.clear();
        mOnes = 0;
        if (isSparse) {
            for (uint32_t i = 0; i < nBits; ++i) {
                if (intArray[i] > 0 != _zero) {
                    indices.insert(i);
                }
            }
            mOnes = indices.size();
            return;
        }
        bitArray = new uint8_t[nBytes];
        for (uint32_t i = 0; i < nBytes; ++i) {
            bitArray[i] = 0;
            for (uint32_t j = 0; j < 8 && i*8+j < nBits; ++j) {
                if (intArray[i*8+j] > 0 != _zero) {
                    bitArray[i] |= (uint8_t) 1 << (7-j);
                    mOnes++;
                }
            }
        }
    }
};

int32_t binary_intersect_mix(const BinaryVariant& lhs, const BinaryVariant& rhs);

int32_t binary_intersect(const BinaryVariant& lhs, const BinaryVariant& rhs);

void parseBinaryFile(const std::string& filename, std::vector<BinaryVariant>& mixStorage);

#endif
