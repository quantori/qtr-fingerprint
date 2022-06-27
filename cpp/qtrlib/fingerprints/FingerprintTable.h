#pragma once

#include "Fingerprint.h"

#include <vector>

namespace qtr {

template<size_t fingerprintSizeInBytes>
using FingerprintTable = std::vector<Fingerprint<fingerprintSizeInBytes>>;

template<size_t fingerprintSizeInBytes>
FingerprintTable<fingerprintSizeInBytes> buildFingerprintTableFromSDFFile(const std::string &sdfFile);

using IndigoFingerprintTable = FingerprintTable<IndigoFingerprint::sizeInBytes>;

} // namespace qtr

#include "FingerprintTable.hpp"