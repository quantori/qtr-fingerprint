#pragma once

#include <vector>
#include <map>

#include "FingerprintTable.h"

namespace qtr {

    std::map<size_t, std::vector<bool>>
    fingerprintsToColumns(const IndigoFingerprintTable &fingerprints, const std::vector<size_t> &columnIds = {});

    std::vector<size_t> sortIndexesByValues(const std::map<size_t, double> &values);

} // namespace qtr
