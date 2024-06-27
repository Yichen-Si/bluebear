#include "binary_variant_compress.h"

int32_t binary_intersect_mix(const BinaryVariant& lhs, const BinaryVariant& rhs) {
    if (!(lhs.isSparse && !rhs.isSparse)) {
        return -1;
    }
    int32_t count = 0;
    for (auto const& idx : lhs.indices) {
        count += (rhs.bitArray[idx >> 3] & (uint8_t) 1 << (7 - idx % 8)) != 0;
    }
    return count;
}

int32_t binary_intersect(const BinaryVariant& lhs, const BinaryVariant& rhs) {
    if (lhs.isSparse && rhs.isSparse) {
        std::vector<uint32_t> intersect;
        auto it = std::set_intersection(lhs.indices.begin(), lhs.indices.end(), rhs.indices.begin(), rhs.indices.end(), std::back_inserter(intersect));
        return intersect.size();
    }
    if (!lhs.isSparse && !rhs.isSparse) {
        return byte_pair_op.sum_and(lhs.bitArray, rhs.bitArray, lhs.nBytes);
    }
    if (rhs.isSparse) {
        return binary_intersect_mix(rhs, lhs);
    }
    return binary_intersect_mix(lhs, rhs);
}

void parseBinaryFile(const std::string& filename, std::vector<BinaryVariant>& mixStorage) {

    std::ifstream infile(filename, std::ios::binary);
    if (!infile) {
        error("Could not open file: %s", filename.c_str());
    }

    int32_t nValues, nalleles, sparse_thres;
    infile.read(reinterpret_cast<char*>(&nValues), sizeof(int32_t));
    infile.read(reinterpret_cast<char*>(&nalleles), sizeof(int32_t));
    infile.read(reinterpret_cast<char*>(&sparse_thres), sizeof(int32_t));
    mixStorage.clear();
    mixStorage.reserve(nValues);

    char delim[2] = {'>', '<'};
    char sparse_indicator[2] = {'0', '1'};
    for (int32_t i = 0; i < nValues; ++i) {
        char delim_read;
        infile.read(&delim_read, 1);
        if (delim_read != delim[0] && delim_read != delim[1]) {
            error("Error in parsing the binary file. (Invalid delimiter %c)", delim_read);
        }
        bool flipped = delim_read == '<';
        int32_t minor_ac, nVariants, npass;
        infile.read(reinterpret_cast<char*>(&minor_ac), sizeof(int32_t));
        infile.read(reinterpret_cast<char*>(&nVariants), sizeof(int32_t));
        infile.read(reinterpret_cast<char*>(&npass), sizeof(int32_t));
        char sparse_read;
        infile.read(&sparse_read, 1);
        if (sparse_read != sparse_indicator[0] &&
            sparse_read != sparse_indicator[1]) {
            error("Error in parsing the binary file. (Invalid sparse indicator %c)", sparse_read);
        }
        mixStorage.emplace_back(sparse_read == '1', nalleles, minor_ac, flipped);
        BinaryVariant& bv = mixStorage.back();
        if (bv.isSparse) {
            for (int j = 0; j < minor_ac; ++j) {
                uint32_t index;
                infile.read(reinterpret_cast<char*>(&index), sizeof(uint32_t));
                bv.indices.insert(index);
            }
        } else {
            bv.bitArray = new uint8_t[bv.nBytes];
            infile.read(reinterpret_cast<char*>(bv.bitArray), bv.nBytes);
        }
    }
    infile.close();
}
