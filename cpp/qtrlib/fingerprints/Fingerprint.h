#pragma once

#include <cstdint>
#include <cstring>
#include "IndigoMolecule.h"
#include "Utils.h"
#include "indigo.h"

template<std::size_t FINGERPRINT_SIZE>
class FingerprintTest;

/**
 * @tparam FINGERPRINT_SIZE fingerprint size in bits
 */
template<std::size_t FINGERPRINT_SIZE>
class __attribute__ ((__packed__)) Fingerprint {
public:
    typedef uint64_t BlockType;

    /**
     * Create fingerprint, data is full with zero
     */
    Fingerprint() {
        memset(_data, 0, countOfBlocks * sizeof(BlockType));
    }

    Fingerprint(indigo_cpp::IndigoMolecule &mol) {
        mol.aromatize();
        int fingerprint = indigoFingerprint(mol.id(), "sub");
        const char *fp = indigoToString(fingerprint);
        buildFromIndigoFingerprint(fp);

    }

    ~Fingerprint() = default;

    Fingerprint(const Fingerprint &other) = default;

    [[nodiscard]] bool get(std::size_t index) const {
        return getBit(_data[index / blockTypeBits],
                      index % blockTypeBits);
    }

private:
    const inline static std::size_t blockTypeBits = fromBytesToBits(sizeof(BlockType));
    const inline static std::size_t countOfBlocks = divideIntegersCeil(FINGERPRINT_SIZE, blockTypeBits);
    BlockType _data[countOfBlocks];

    void buildFromIndigoFingerprint(const char *fp) {
        const int CNT_BIT_HEX = 4;
        for (int i = 0; fp[i] != '\0'; ++i) {
            BlockType intNumber = (fp[i] >= 'a') ? (fp[i] - 'a' + 10) : (fp[i] - '0');
            _data[(CNT_BIT_HEX * i) / blockTypeBits] += intNumber << (blockTypeBits - CNT_BIT_HEX * (i + 1));
        }
    }

    friend FingerprintTest<FINGERPRINT_SIZE>;
};