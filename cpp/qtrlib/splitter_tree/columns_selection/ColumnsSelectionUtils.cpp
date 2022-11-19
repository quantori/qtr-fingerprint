#include "ColumnsSelectionUtils.h"

#include <algorithm>
#include <map>

namespace qtr {

    namespace {
        std::map<size_t, std::vector<bool>> fingerprintToFullColumns(const IndigoFingerprintTable &fingerprints) {
            size_t columnsNumber = fingerprints[0].size();
            std::map<size_t, std::vector<bool>> columns;
            for (auto &fp: fingerprints) {
                for (size_t col = 0; col < columnsNumber; col++) {
                    columns[col].emplace_back(fp[col]);
                }
            }
            return columns;
        }

        std::map<size_t, std::vector<bool>>
        fingerprintToColumnsSubset(const IndigoFingerprintTable &fingerprints, const std::vector<size_t> &columnIds) {
            std::map<size_t, std::vector<bool>> columns;
            for (auto &fp: fingerprints) {
                for (size_t col: columnIds) {
                    columns[col].emplace_back(fp[col]);
                }
            }
            return columns;
        }
    }

    std::map<size_t, std::vector<bool>>
    fingerprintsToColumns(const IndigoFingerprintTable &fingerprints, const std::vector<size_t> &columnIds) {
        if (columnIds.empty()) {
            return fingerprintToFullColumns(fingerprints);
        } else {
            return fingerprintToColumnsSubset(fingerprints, columnIds);
        }
    }

    std::vector<size_t> sortIndexesByValues(const std::map<size_t, double> &values) {
        std::vector<size_t> indexes;
        for (auto &[index, _]: values) {
            indexes.emplace_back(index);
        }
        std::sort(indexes.begin(), indexes.end(), [&values](size_t a, size_t b) {
            return values.find(a)->second < values.find(b)->second;
        });
        return indexes;
    }
} // namespace qtr