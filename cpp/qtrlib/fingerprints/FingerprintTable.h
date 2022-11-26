#pragma once

#include "Fingerprint.h"

#include <vector>

namespace qtr {

    template<size_t fingerprintSizeInBits>
    class FingerprintTable : public std::vector<Fingerprint<fingerprintSizeInBits>> {
    public:
        FingerprintTable() = default;

        explicit FingerprintTable(const std::string &sdfFile);
    };

    using IndigoFingerprintTable = FingerprintTable<IndigoFingerprint::size()>;

} // namespace qtr
