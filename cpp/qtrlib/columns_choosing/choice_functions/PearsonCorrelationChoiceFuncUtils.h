#pragma once

#include <vector>
#include "FingerprintTable.h"

namespace qtr {

    constexpr static const double pearsonCorrelationInfinityValue = 2.0;

    double findPearsonCorrelation(const std::vector<bool> &x, const std::vector<bool> &y);

    bool isConstColumn(const std::vector<bool> &x);

    std::vector<double> findMaxAbsPearsonCorrelation(const std::vector<std::vector<bool>> &columns);

    qtr::IndigoFingerprintTable chooseSubset(const IndigoFingerprintTable& fingerprints, size_t subsetSize);

} // namespace qtr
