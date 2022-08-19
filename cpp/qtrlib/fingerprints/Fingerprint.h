#pragma once

#include <bitset>
#include <climits>
#include <cstddef>
#include <vector>

#include "Utils.h"

class QtrIndigoFingerprint;

namespace qtr {

template<size_t fingerprintSizeInBytes>
class Fingerprint : public std::bitset<CHAR_BIT*fingerprintSizeInBytes>
{
public:
    static constexpr size_t sizeInBytes = fingerprintSizeInBytes;
    static constexpr size_t sizeInBits = fromBytesToBits(fingerprintSizeInBytes);

    Fingerprint() = default;

    explicit Fingerprint(const QtrIndigoFingerprint &f);

    /**
     * Builds fingerprint from HEX string(for example - "a4210f")
     * @param s - string with fingerprint, strlen(s) should be equals to fingerprintSizeInBytes * 2
     */
    explicit Fingerprint(const std::string &s) {
        for (size_t i = 0; i < s.size(); ++i) {
            size_t j = i * 4ull;
            int currentSym = chexToInt(s[i]);
            this->operator[](j + 3) = currentSym & 1;
            this->operator[](j + 2) = currentSym & 2;
            this->operator[](j + 1) = currentSym & 4;
            this->operator[](j + 0) = currentSym & 8;
        }
    }

    void setBytes(const std::vector<std::byte> &bytes) {
        this->reset();
        for(size_t i = 0; i < bytes.size(); i++)
            for(size_t j = 0; j < CHAR_BIT; j++)
                if (bool((bytes[i] >> j) & std::byte(1))) 
                    this->set(i*CHAR_BIT + j);
    }

    std::vector<std::byte> getBytes() const {
        std::vector<std::byte> result(fingerprintSizeInBytes, std::byte(0));

        for(size_t i = 0; i < result.size(); i++)
            for(size_t j = 0; j < CHAR_BIT; j++)
                if (this->test(i*CHAR_BIT + j)) 
                    result[i] |= (std::byte(1) << j);

        return result;
    }

    void saveBytes(std::ostream& out) const {
        auto toWrite = getBytes();
        out.write((char*)toWrite.data(), toWrite.size());
    }

    void readFrom(std::istream& in) {
        for (uint64_t i = 0, j = 0;
             i < sizeInBytes; ++i, j += BIT_IN_BYTE) {
            auto curr = in.get();
            for (uint64_t k = 0; k < BIT_IN_BYTE; ++k)
                this->operator[](j + k) = (curr & (1ull << k));
        }
    }
};

using IndigoFingerprint = Fingerprint<323>; // todo set size of fingerprint by number of bits, not bytes
using FullIndigoFingerprint = Fingerprint<467>;

} // namespace qtr