#pragma once

#include <cstdint>
#include <cstring>

#include "Utils.h"

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

    ~Fingerprint() = default;

    Fingerprint(const Fingerprint &other) = default;

    [[nodiscard]] bool get(std::size_t index) const {
        return getBit(_data[index / blockTypeBits],
                      index % blockTypeBits);
    }

private:
    const static std::size_t blockTypeBits = fromBytesToBits(sizeof(BlockType));
    const static std::size_t countOfBlocks = divideIntegersCeil(FINGERPRINT_SIZE, blockTypeBits);
    BlockType _data[countOfBlocks];
};