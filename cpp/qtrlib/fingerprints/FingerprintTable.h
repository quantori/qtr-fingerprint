#pragma once

#include "Fingerprint.h"

#include <vector>

namespace qtr {

template<size_t fingerprintSizeInBytes>
class FingerprintTable : public std::vector<Fingerprint<fingerprintSizeInBytes>> {
public:
    FingerprintTable() = default;

    explicit FingerprintTable(const std::string &sdfFile);
};

using IndigoFingerprintTable = FingerprintTable<IndigoFingerprint::sizeInBytes>;

} // namespace qtr
