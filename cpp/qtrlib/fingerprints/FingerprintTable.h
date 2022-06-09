#pragma once

#include "Fingerprint.h"

#include <vector>

namespace qtr {

template<size_t fingerprintSizeInBytes>
using FingerprintTable = std::vector<Fingerprint<fingerprintSizeInBytes>>;

using IndigoFingerprintTable = FingerprintTable<IndigoFingerprint::sizeInBytes>;

} // namespace qtr