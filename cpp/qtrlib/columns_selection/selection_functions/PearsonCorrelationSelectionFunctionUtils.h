#pragma once

#include <vector>
#include <map>
#include "FingerprintTable.h"
#include "columns_io/ColumnsIOConsts.h"

namespace qtr {

    constexpr static const double pearsonCorrelationInfinityValue = 2.0;

    double findPearsonCorrelation(const std::vector<bool> &x, const std::vector<bool> &y);

    bool isConstColumn(const std::vector<bool> &x);

    std::map<size_t, double> findMaxAbsPearsonCorrelation(const std::map<size_t, std::vector<bool>> &columns);

    qtr::IndigoFingerprintTable selectFingerprintSubset(const IndigoFingerprintTable &fingerprints, size_t subsetSize);

} // namespace qtr
