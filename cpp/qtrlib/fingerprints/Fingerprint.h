#pragma once

#include <cstdint>
#include <cstring>
#include "IndigoMolecule.h"
#include "Utils.h"
#include "indigo.h"

const inline int CNT_BIT_HEX = 4; // cnt bits in one hex number

template<std::size_t FINGERPRINT_SIZE>
class FingerprintTest;

class FingerprintTableTest;

/**
 * @tparam FINGERPRINT_SIZE fingerprint size in bits
 */
template<std::size_t FINGERPRINT_SIZE>
class __attribute__ ((__packed__)) Fingerprint {
public:
    typedef uint64_t BlockType;
    const inline static std::size_t blockTypeBits = fromBytesToBits(sizeof(BlockType));
    const inline static std::size_t countOfBlocks = divideIntegersCeil(FINGERPRINT_SIZE, blockTypeBits);

    /**
     * Create fingerprint, data is full with zero
     */
    Fingerprint() {
        memset(_data, 0, countOfBlocks * sizeof(BlockType));
    }

    Fingerprint(indigo_cpp::IndigoMolecule &mol);

    ~Fingerprint() = default;

    Fingerprint(const Fingerprint &other) = default;

    [[nodiscard]] bool get(std::size_t index) const {
        return getBit(_data[index / blockTypeBits],
                      (blockTypeBits - 1 - index) % blockTypeBits);
    }

    int lenInBinary(const int &indexHex) {
        return CNT_BIT_HEX * indexHex;
    }

    int getBlockNumber(const int &indexHex) {
        return lenInBinary(indexHex) / blockTypeBits;
    }

private:
    BlockType _data[countOfBlocks];

    void buildFromIndigoFingerprint(const char *fp);

    friend FingerprintTest<FINGERPRINT_SIZE>;
    friend FingerprintTableTest;
};

#include "Fingerprints.hpp"
