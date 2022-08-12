#include "PearsonCorrelationSelectionFunctionUtils.h"
#include "columns_io/ColumnsIOConsts.h"

#include <cassert>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <random>

namespace qtr {

    namespace {
        std::mt19937 colSelectionRandom(123);
    }

    bool isConstColumn(const std::vector<bool> &x) {
        return std::all_of(x.begin(), x.end(), [val = x[0]](bool a) {
            return a == val;
        });
    }

    double findPearsonCorrelation(const std::vector<bool> &x, const std::vector<bool> &y) {
        assert(x.size() == y.size());
        int64_t n = (int32_t) x.size();
        int64_t c[2][2] = {0, 0, 0, 0};
        int64_t xSum = 0, ySum = 0;
        for (size_t i = 0; i < x.size(); i++) {
            c[x[i]][y[i]]++;
            xSum += x[i];
            ySum += y[i];
        }

        if (xSum == x.size() || xSum == 0 || ySum == y.size() || ySum == 0)
            throw std::invalid_argument("x and y can not be constant columns");

        int64_t numerator = 0;
        for (int64_t i = 0; i <= 1; i++) {
            for (int64_t j = 0; j <= 1; j++) {
                numerator += c[i][j] * (n * i - xSum) * (n * j - ySum);
            }
        }
        int64_t denominator = xSum * (n - xSum) * ySum * (n - ySum);
        return double(numerator) / (double) n / sqrt((double) denominator);
    }

    std::vector<double> findMaxAbsPearsonCorrelation(const std::vector<std::vector<bool>> &columns,
                                                     const std::vector<size_t> &columnsIndexes = {}) {
        std::vector<double> maxCorrelation(columns.size(), pearsonCorrelationInfinityValue);
        const std::vector<size_t> *indexes = &columnsIndexes;
        std::vector<size_t> fullIndexes;
        if (columnsIndexes.empty()) {
            fullIndexes.resize(columns.size());
            std::iota(fullIndexes.begin(), fullIndexes.end(), 0);
            indexes = &fullIndexes;
        }
        for (size_t i = 0; i < indexes->size(); i++) {
            size_t x = (*indexes)[i];
            if (isConstColumn(columns[x])) {
                maxCorrelation[x] = pearsonCorrelationInfinityValue;
                continue;
            }
            for (size_t j = i + 1; j < indexes->size(); j++) {
                size_t y = (*indexes)[j];
                if (isConstColumn(columns[y]))
                    continue;
                double corr = std::abs(findPearsonCorrelation(columns[x], columns[y]));
                maxCorrelation[x] = std::max(maxCorrelation[x], corr);
                maxCorrelation[y] = std::max(maxCorrelation[y], corr);
            }
        }
        return maxCorrelation;
    }

    qtr::IndigoFingerprintTable selectFingerprintSubset(const IndigoFingerprintTable &fingerprints, size_t subsetSize) {
        std::vector<size_t> indexes(fingerprints.size());
        std::iota(indexes.begin(), indexes.end(), 0);
        std::shuffle(indexes.begin(), indexes.end(), colSelectionRandom);
        indexes.resize(std::min(subsetSize, indexes.size()));
        IndigoFingerprintTable subset;
        for (size_t index: indexes) {
            subset.emplace_back(fingerprints[index]);
        }
        return subset;
    }

} // namespace qtr
