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

#include <DataStructs/ExplicitBitVect.h>

#include "Utils.h"
#include "Bitset.h"
#include "QtrIndigoFingerprint.h"


namespace qtr {

    class Fingerprint : public qtr::Bitset<> {
    public:
        using Bitset<>::sizeInBytes;
        using Bitset<>::size;

        explicit Fingerprint(size_t size);

        Fingerprint(const Fingerprint &other) = default;

        Fingerprint() = default;

        /**
         * Builds fingerprint from HEX string(for example - "a4210f")
         * @param s - string with fingerprint, strlen(s) should be equals to fingerprintSizeInBytes * 2
         */
        explicit Fingerprint(const std::string &s);

        /**
         * Builds fingerprint from byte array
         * @param arr
         */
        explicit Fingerprint(const indigo::Array<byte> &arr);

        explicit Fingerprint(const QtrIndigoFingerprint &f);

        explicit Fingerprint(std::unique_ptr<ExplicitBitVect> f);

        void setBytes(const std::vector<std::byte> &bytes);

        [[nodiscard]] std::vector<std::byte> getBytes() const;

    private:
        void addSymbol(size_t blockIndex, size_t blockSize, size_t symbol);
    };

    constexpr size_t IndigoFingerprintSize = 2584;
    constexpr size_t FullIndigoFingerprintSize = 3736;

    static constexpr std::pair<uint64_t, uint64_t> FullIndigoFingerprintEmptyRange = {1624, 2776};

    static_assert(FullIndigoFingerprintSize - IndigoFingerprintSize ==
                  FullIndigoFingerprintEmptyRange.second - FullIndigoFingerprintEmptyRange.first);

    Fingerprint cutFullIndigoFingerprint(const Fingerprint &fullFingerprint);

    Fingerprint indigoFingerprintFromSmiles(const std::string &smiles);

    Fingerprint rdkitFingerprintFromSmiles(const std::string &smiles);

} // namespace qtr