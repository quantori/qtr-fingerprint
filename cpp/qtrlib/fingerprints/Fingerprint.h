#pragma once

#include <bitset>
#include <climits>
#include <cstddef>
#include <vector>
#include <cassert>
#include <string>

#include "indigo.h"

#include "IndigoMolecule.h"
#include "IndigoSession.h"
#include "IndigoWriteBuffer.h"

#include "Utils.h"
#include "Bitset.h"
#include "QtrIndigoFingerprint.h"

namespace qtr {

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
            assert(s.size() == sizeInBytes * 2);
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
            assert(bytes.size() == sizeInBytes);
            for (size_t i = 0; i < bytes.size(); i++)
                for (size_t j = 0; j < CHAR_BIT; j++)
                    this->operator[](i * CHAR_BIT + j) = bool((bytes[i] >> j) & std::byte(1));
        }

        [[nodiscard]] std::vector<std::byte> getBytes() const {
            std::vector<std::byte> result(sizeInBytes, std::byte(0));
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

    static constexpr std::pair<uint64_t, uint64_t> FullIndigoFingerprintEmptyRange = {1624, 2776};

    static_assert(FullIndigoFingerprint::size() - IndigoFingerprint::size() ==
                  FullIndigoFingerprintEmptyRange.second - FullIndigoFingerprintEmptyRange.first);

    inline IndigoFingerprint cutFullIndigoFingerprint(const FullIndigoFingerprint &fullFingerprint) {
        IndigoFingerprint fingerprint;
        for (size_t i = 0; i < FullIndigoFingerprint::size(); i++) {
            if (i < FullIndigoFingerprintEmptyRange.first) {
                fingerprint[i] = fullFingerprint[i];
            } else if (i >= FullIndigoFingerprintEmptyRange.second) {
                fingerprint[i - FullIndigoFingerprintEmptyRange.second +
                            FullIndigoFingerprintEmptyRange.first] = fullFingerprint[i];
            } else {
                assert(!fullFingerprint[i]);
            }
        }
        return fingerprint;
    }

    inline IndigoFingerprint indigoFingerprintFromSmiles(const std::string &smiles) {
        auto indigoSessionPtr = indigo_cpp::IndigoSession::create();
        auto mol = indigoSessionPtr->loadMolecule(smiles);
        mol.aromatize();
        int fingerprint = indigoFingerprint(mol.id(), "sub");
        FullIndigoFingerprint fullFingerprints(indigoToString(fingerprint));
        IndigoFingerprint cutFingerprint = cutFullIndigoFingerprint(fullFingerprints);
        return cutFingerprint;
    }
} // namespace qtr