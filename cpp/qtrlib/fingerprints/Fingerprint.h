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

#include "base_cpp/array.h"
#include "base_cpp/scanner.h"
#include "src/bingo_object.h"
#include "src/indigo_internal.h"
#include "molecule/smiles_loader.h"
#include "molecule/molecule_fingerprint.h"

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
                addSymbol(i, 4, (size_t) chexToInt(s[i]));
            }
        }

        /**
         * Builds fingerprint from byte array
         * @param arr
         */

        explicit Fingerprint(const indigo::Array<byte> &arr) {
            assert(arr.size() == sizeInBytes);
            for (size_t i = 0; i < arr.size(); i++) {
                addSymbol(i, 8, arr[i]);
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

    private:
        void addSymbol(size_t blockIndex, size_t blockSize, size_t symbol) {
            for (size_t i = 0; i < blockSize; i++) {
                (*this)[blockIndex * blockSize + blockSize - i - 1] = symbol & (1u << i);
            }
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

    const static Indigo indigoInstance;

    inline IndigoFingerprint indigoFingerprintFromSmiles(const std::string &smiles) {
        indigo::BufferScanner scanner(smiles.c_str(), smiles.size(), false);
        indigo::SmilesLoader loader(scanner);
        indigo::Molecule molecule;
        loader.loadMolecule(molecule);
        molecule.aromatize(indigo::AromaticityOptions());
        bingo::IndexMolecule indexMolecule(molecule, indigo::AromaticityOptions());
        indigo::Array<byte> subFingerprint;
        indigo::MoleculeFingerprintBuilder fingerprintBuilder(molecule, indigoInstance.fp_params);
        fingerprintBuilder.parseFingerprintType("sub", false);
        fingerprintBuilder.process();
        subFingerprint.copy(fingerprintBuilder.get(), indigoInstance.fp_params.fingerprintSize());
//        indexMolecule.buildFingerprint(indigoInstance.fp_params, &subFingerprint, nullptr);
        FullIndigoFingerprint fullIndigoFingerprint(subFingerprint);
        IndigoFingerprint  cutFingerprint = cutFullIndigoFingerprint(fullIndigoFingerprint);
        return cutFingerprint;

        // TODO: implement indigoFingerprint Building using indigo core

//        auto indigoSessionPtr = indigo_cpp::IndigoSession::create();
//        auto mol = indigoSessionPtr->loadMolecule(smiles);
//        mol.aromatize();
//        int fingerprint = indigoFingerprint(mol.id(), "sub");
//        FullIndigoFingerprint fullFingerprints(indigoToString(fingerprint));
//        IndigoFingerprint cutFingerprint = cutFullIndigoFingerprint(fullFingerprints);
//        return cutFingerprint;
    }
} // namespace qtr