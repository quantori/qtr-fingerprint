#pragma once

#include <vector>

namespace qtr {

    double findPearsonCorrelation(const std::vector<bool> &x, const std::vector<bool> &y);

    std::vector<double> findMaxAbsPearsonCorrelation(const std::vector<std::vector<bool>> &columns);

} // namespace qtr
