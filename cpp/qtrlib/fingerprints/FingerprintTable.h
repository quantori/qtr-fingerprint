#pragma once

#include "Fingerprint.h"

#include <vector>

template<size_t fingerprintSizeInBytes>
using FingerprintTable = std::vector<Fingerprint<fingerprintSizeInBytes>>;

using FingerprintTableForIndigo = FingerprintTable<FingerprintForIndigo::sizeInBytes>;