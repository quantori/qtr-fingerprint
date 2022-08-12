#pragma once

#include <vector>
#include "FingerprintTable.h"
#include "columns_io/ColumnsIOConsts.h"

namespace qtr {

    constexpr static const double pearsonCorrelationInfinityValue = 2.0;

    double findPearsonCorrelation(const std::vector<bool> &x, const std::vector<bool> &y);

    bool isConstColumn(const std::vector<bool> &x);

    std::vector<double> findMaxAbsPearsonCorrelation(const std::vector<std::vector<bool>> &columns,
                                                     const std::vector<size_t> &columnsIndexes);

    qtr::IndigoFingerprintTable selectFingerprintSubset(const IndigoFingerprintTable &fingerprints, size_t subsetSize);

} // namespace qtr
