#include <cassert>
#include <cmath>
#include <numeric>

#include "PearsonCorrelationChoiceFunc.h"
#include "FingerprintTable.h"

namespace qtr {

    static double
    findCorrelation(const std::vector<bool> &x, const std::vector<bool> &y) {
        assert(x.size() == y.size());
        size_t c[2][2] = {0, 0, 0, 0};
        uint32_t xSum = 0, ySum = 0;
        for (size_t i = 0; i < x.size(); i++) {
            c[x[i]][y[i]]++;
            xSum += x[i];
            ySum += y[i];
        }
        double xMean = (double) xSum / (double) x.size();
        double yMean = (double) ySum / (double) y.size();
        double numerator = 0;
        for (size_t i = 0; i <= 1; i++) {
            for (size_t j = 0; j <= 1; j++) {
                numerator += (double) c[i][j] * ((double) i - xMean) * ((double) j - yMean);
            }
        }
        double denominator = (xSum * (1 - xMean) * (1 - xMean) + (1 - xSum) * xMean * xMean) *
                             (ySum * (1 - yMean) * (1 - yMean) + (1 - yMean) * yMean * yMean);
        return numerator / sqrt(denominator);
    }

    static std::vector<double>
    findMaxCorrelation(const std::vector<std::vector<bool>> &columns) {
        std::vector<double> maxCorrelation(columns.size());
        for (size_t col1 = 0; col1 < columns.size(); col1++) {
            for (size_t col2 = col1 + 1; col2 < columns.size(); col2++) {
                double corr = std::abs(findCorrelation(columns[col1], columns[col2]));
                maxCorrelation[col1] = std::max(maxCorrelation[col1], corr);
                maxCorrelation[col2] = std::max(maxCorrelation[col2], corr);
            }
        }
        return maxCorrelation;
    }

    static std::vector<std::vector<bool>> fingerprintsToColumns(const IndigoFingerprintTable &fingerprints) {
        size_t columnsNumber = fingerprints[0].size();
        std::vector<std::vector<bool>> columns(fingerprints.size());
        for (auto &fp: fingerprints) {
            for (size_t col = 0; col < columnsNumber; col++) {
                columns[col].emplace_back(fp[col]);
            }
        }
        return columns;
    }

    static std::vector<int> findOrder(const std::vector<double> &maxCorrelation) {
        std::vector<int> order(maxCorrelation.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&maxCorrelation](int a, int b) {
            return maxCorrelation[a] < maxCorrelation[b];
        });
        return order;
    }

    choice_result_t PearsonCorrelationChoiceFunc::operator()(const IndigoFingerprintTable &fingerprints) const {
        assert(!fingerprints.empty());
        auto columns = fingerprintsToColumns(fingerprints);
        auto maxCorrelation = findMaxCorrelation(columns);
        return findOrder(maxCorrelation);
    }

} // namespace qtr