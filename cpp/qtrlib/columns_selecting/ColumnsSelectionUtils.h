#pragma once

#include <vector>

#include "FingerprintTable.h"

namespace qtr {

    std::vector<std::vector<bool>> fingerprintsToColumns(const IndigoFingerprintTable &fingerprints);

    std::vector<size_t> sortIndexesByValues(const std::vector<double> &values);

} // namespace qtr
