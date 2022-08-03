#include "ColumnsChoiceUtils.h"

#include <numeric>
#include <algorithm>

namespace qtr {
    std::vector<std::vector<bool>> fingerprintsToColumns(const IndigoFingerprintTable &fingerprints) {
        size_t columnsNumber = fingerprints[0].size();
        std::vector<std::vector<bool>> columns(columnsNumber);
        for (auto &fp: fingerprints) {
            for (size_t col = 0; col < columnsNumber; col++) {
                columns[col].emplace_back(fp[col]);
            }
        }
        return columns;
    }

    std::vector<int> sortIndexesByValues(const std::vector<double> &values) {
        std::vector<int> indexes(values.size());
        std::iota(indexes.begin(), indexes.end(), 0);
        std::sort(indexes.begin(), indexes.end(), [&values](int a, int b) {
            return values[a] < values[b];
        });
        return indexes;
    }
} // namespace qtr