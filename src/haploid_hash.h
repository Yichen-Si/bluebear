#ifndef __HAPLOID_HASH_H
#define __HAPLOID_HASH_H

#include <set>
#include <tuple>
#include <iostream>
#include <vector>
#include <functional>
#include <unordered_map>
#include <random>
#include <algorithm>
#include <cstring>
#include "byte_operation.h"
#include "binary_variant_compress.h"

/**
 * Hash a binary array or a set of indices to a 64-bit integer
 * Assuming the array represent an often sparse binary variant
 * Currently - use the indices of the first positions of consecutive 1s
*/
class LshBinary {
public:
    uint32_t nBits, nBytes;
    int32_t nValues;
    int32_t nMarkers; // Number of change positions to consider
    int32_t nQuery, nCompare, nCollision;
    std::vector<int32_t> nVariants; // Number of input elements with exactly the same hash key
    uint64_t magicConst64 = 0x9e3779b97f4a7c15;
    uint64_t magicConst32 = 0x9e3779b9;

    LshBinary(uint32_t _n, int32_t _m) : nBits(_n), nMarkers(_m) {
        nBytes = (nBits + 7) / 8;
        nValues = 0;
        nQuery = 0;
        nCompare = 0;
        nCollision = 0;
    }

    // Not sure if this is the best way
    // We do not spread the influence of each additional value
    // Assuming most of the times we only need to hash a few 16~32 bit values
    void combine_hash(uint64_t& key, uint64_t value) {
        key ^= value + magicConst64;
        key = (key << 16) | (key >> 48);
    }

    // Hash based on the indices of 1s (or 0s if hashZeros) in the bitArray
    uint64_t hash_bit(const uint8_t* bitArray, bool hashZeros = false, int32_t maxIdx = 0) {
        uint64_t key = 0;
        maxIdx = (maxIdx > 0) ? std::min(maxIdx, nMarkers) : nMarkers;
        int32_t count = 0;
        bool prev_bit = 0;
        for (uint32_t i = 0; i < nBytes; ++i) {
            for (uint32_t j = 0; j < 8 && i*8+j < nBits; ++j) {
                bool current_bit = (bitArray[i] & (1 << (7 - j))) > 0 != hashZeros;
                if (current_bit && (!prev_bit)) {
                    combine_hash(key, static_cast<uint64_t>(i*8+j));
                    count++;
                    if (count > maxIdx) {
                        return key;
                    }
                }
                prev_bit = current_bit;
            }
        }
        return key;
    }

    template <typename T>
    uint64_t hash_int(const T* intArray, bool hashZeros = 0, int32_t maxIdx = 0) {
        uint64_t key = 0;
        maxIdx = (maxIdx > 0) ? std::min(maxIdx, nMarkers) : nMarkers;
        int32_t count = 0;
        bool prev_bit = 0;
        for (uint32_t i = 0; i < nBits; ++i) {
            bool current_bit = (intArray[i] > 0) != hashZeros;
            if (current_bit && (!prev_bit)) {
                combine_hash(key, static_cast<uint64_t>(i));
                count++;
                if (count > maxIdx) {
                    return key;
                }
            }
            prev_bit = current_bit;
        }
        return key;
    }

    uint64_t hash_set(const std::set<uint32_t>& indices, int32_t maxIdx = 0) {
        uint64_t key = 0;
        maxIdx = (maxIdx > 0) ? std::min(maxIdx, nMarkers) : nMarkers;
        int32_t count = 0;
        for (auto index : indices) { // rely on the fact that the set is sorted
            combine_hash(key, static_cast<uint64_t>(index));
            count++;
            if (count > maxIdx) {
                return key;
            }
        }
        return key;
    }
};

class LshBinaryMix : public LshBinary {
public:
    std::unordered_map<uint64_t, std::vector<size_t>> table; // hash key -> indices
    std::vector<BinaryVariant> mixStorage;
    bool oriented;   // If oriented, 01000 and 10111 are different
    int32_t sparseThreshold;

    LshBinaryMix(uint32_t _nBits, int32_t _nMarkers, int32_t _sThres, bool _o = true)
    : LshBinary(_nBits, _nMarkers), sparseThreshold(_sThres), oriented(_o) {}

    int32_t find_exact(uint64_t key, BinaryVariant &variant) {
        auto ptr = table.find(key);
        if (ptr == table.end()) {
            return -1;
        }
        for (auto index : ptr->second) {
            nCompare++;
            if (mixStorage[index].equal(variant)) {
                nCollision++;
                return index;
            }
        }
        return -1;
    }

    int32_t find_exact(uint64_t key, std::set<uint32_t>& indices, bool sign = 0) {
        if(indices.size() > sparseThreshold) { // should never happen
            return -2;
        }
        auto ptr = table.find(key);
        if (ptr == table.end()) {
            return -1;
        }
        for (auto index : ptr->second) {
            if (!mixStorage[index].isSparse) {
                continue;
            }
            if (oriented && mixStorage[index].sign != sign) {
                continue;
            }
            nCompare++;
            if (mixStorage[index].indices == indices) {
                nCollision++;
                return index;
            }
        }
        return -1;
    }

    int32_t find_exact(uint64_t key, const uint8_t* bitArray, bool sign = 0) {
        auto ptr = table.find(key);
        if (ptr == table.end()) {
            return -1;
        }
        for (auto index : ptr->second) {
            if (mixStorage[index].isSparse) {
                continue;
            }
            nCompare++;
            if (oriented && mixStorage[index].sign != sign) {
                bool equal = true;
                for (uint32_t i = 0; i < nBytes - 1; ++i) {
                    if (mixStorage[index].bitArray[i] != ~bitArray[i]) {
                        equal = false;
                        break;
                    }
                }
                if (equal) {
                    uint8_t mask = 255 << (8 - nBits % 8);
                    if ( (mixStorage[index].bitArray[nBytes-1] & mask) == ((~bitArray[nBytes-1]) & mask) ) {
                        nCollision++;
                        return index;
                    }
                }
            } else if (memcmp(mixStorage[index].bitArray, bitArray, nBytes) == 0) {
                nCollision++;
                return index;
            }
        }
        return -1;
    }

    template <typename T>
    int32_t add_int(const T* intArray, bool isSparse, bool zero, int32_t maxIdx = 0) {
        nQuery++;
        bool sign = oriented && zero;
        BinaryVariant variant(isSparse, nBits, 0, sign);
        variant.update_int<T>(intArray, zero);
        uint64_t key = hash_int<T>(intArray, zero, maxIdx);
        int32_t index = find_exact(key, variant);
        if (index >= 0) {
            nVariants[index]++;
            return index;
        }
        mixStorage.push_back(std::move(variant));
        table.emplace(key, std::vector<size_t>{(size_t) nValues});
        nVariants.push_back(1);
        nValues++;
        return - 1;
    }

    int32_t add_set(std::set<uint32_t>& indices, bool _o = 0) {
        nQuery++;
        bool sign = oriented && _o;
        uint64_t key = hash_set(indices);
        int32_t index = find_exact(key, indices, sign);
        if (index >= 0) {
            nVariants[index]++;
            return index;
        }
        mixStorage.emplace_back(indices, nBits, sign);
        table[key].push_back(nValues);
        nVariants.push_back(1);
        nValues++;
        return - 1;
    }

};

// Need to test this
class LshBinarySparse : public LshBinary {
public:
    std::unordered_map<uint64_t, std::vector<size_t>> table;
    std::vector<std::set<uint32_t> > storage;

    LshBinarySparse(int32_t nBits, int32_t nMarkers) : LshBinary(nBits, nMarkers) {}

    int32_t find_exact(uint64_t key, std::set<uint32_t>& values) {
        if (table.find(key) == table.end()) {
            return -1;
        }
        for (auto index : table[key]) {
            nCompare++;
            if (storage[index] == values) {
                nCollision++;
                return index;
            }
        }
        return -1;
    }

    int32_t add_set(std::set<uint32_t>& values, int32_t maxIdx = 0) {
        nQuery++;
        uint64_t key = hash_set(values, maxIdx);
        int32_t index = find_exact(key, values);
        if (index >= 0) {
            nVariants[index]++;
            return index;
        }
        storage.push_back(std::move(values));
        table[key].push_back(nValues);
        nVariants.push_back(1);
        nValues++;
        return nValues - 1;
    }

    int32_t add_bit(const uint8_t* bitArray, bool hashZeros = false, int32_t maxIdx = 0) {
        // bit to set
        std::set<uint32_t> indices;
        uint8_t zero = hashZeros ? 1 : 0;
        for (uint32_t i = 0; i < nBytes; ++i) {
            for (uint32_t j = 0; j < 8; ++j) {
                if ((bitArray[i] & (1 << (7 - j))) != zero) {
                    indices.insert(i*8+j);
                }
            }
        }
        return add_set(indices, maxIdx);
    }

    template <typename T>
    int32_t add_int(const T* intArray, bool hashZeros = false, int32_t maxIdx = 0) {
        std::set<uint32_t> indices;
        for (uint32_t i = 0; i < nBits; ++i) {
            if ((intArray[i] > 0) != hashZeros) {
                indices.insert(i);
            }
        }
        return add_set(indices, maxIdx);
    }
};

#endif
