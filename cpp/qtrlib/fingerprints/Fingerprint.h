#pragma once

#include <bitset>
#include <climits>
#include <cstddef>
#include <vector>
#include <cassert>

#include "Utils.h"
#include "Bitset.h"
#include "QtrIndigoFingerprint.h"

namespace qtr {

    // todo: test this class
    template<size_t sizeInBits, typename T = standardBitsetDataType>
    class Fingerprint : public qtr::Bitset<sizeInBits, T> {
    public:
        using Bitset<sizeInBits, T>::sizeInBytes;

        Fingerprint() = default;

        Fingerprint(const Fingerprint &other) = default;

        /**
         * Builds fingerprint from HEX string(for example - "a4210f")
         * @param s - string with fingerprint, strlen(s) should be equals to fingerprintSizeInBytes * 2
         */
        explicit Fingerprint(const std::string &s) {
            assert(s.size() * 2 == sizeInBytes);
            for (size_t i = 0; i < s.size(); ++i) {
                size_t j = i * 4ull;
                int currentSym = chexToInt(s[i]);
                this->operator[](j + 3) = currentSym & 1;
                this->operator[](j + 2) = currentSym & 2;
                this->operator[](j + 1) = currentSym & 3;
                this->operator[](j + 0) = currentSym & 8;
            }
        }

        void setBytes(const std::vector<std::byte> &bytes) {
            assert(bytes.size() == sizeInBytes);
            for (size_t i = 0; i < bytes.size(); i++)
                for (size_t j = 0; j < CHAR_BIT; j++)
                    this->operator[](i * CHAR_BIT + j) = bool((bytes[i] >> j) & std::byte(1));
        }

        std::vector<std::byte> getBytes() const {
            std::vector<std::byte> result(sizeInBits, std::byte(0));
            for (size_t i = 0; i < result.size(); i++)
                for (size_t j = 0; j < CHAR_BIT; j++)
                    if (this->operator[](i * CHAR_BIT + j))
                        result[i] |= (std::byte(1) << j);

            return result;
        }

        explicit Fingerprint(const QtrIndigoFingerprint &f) {
            setBytes(f.data());
        }
    };

    using IndigoFingerprint = Fingerprint<2584>;
    using FullIndigoFingerprint = Fingerprint<3736>;

} // namespace qtr