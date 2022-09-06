#include <cmath>
#include <future>

#include "MaxDispersionBitSelector.h"
#include "Fingerprint.h"
#include "fingerprint_table_io/FingerprintTableReader.h"

namespace qtr {

    namespace {
        std::pair<std::vector<size_t>, size_t> findColumnsSumInFile(const std::filesystem::path &filePath) {
            std::vector<size_t> columnsSum(qtr::IndigoFingerprint::size(), 0);
            size_t tableSize = 0;
            for (const auto &[_, fingerprint]: FingerprintTableReader(filePath)) {
                tableSize++;
                for (size_t i = 0; i < IndigoFingerprint::size(); i++)
                    columnsSum[i] += fingerprint[i];
            }
            return {columnsSum, tableSize};
        }

        void mergeResults(std::vector<uint64_t> &currColumnsSum, uint64_t &currTableSize,
                          const std::vector<uint64_t> &fileColumnsSum, uint64_t fileTableSize) {
            currTableSize += fileTableSize;
            for (size_t i = 0; i < IndigoFingerprint::size(); i++) {
                currColumnsSum[i] += fileColumnsSum[i];
            }
        }

        size_t selectBitByColumnsSum(const std::vector<size_t> &columnsSum, size_t tableSize) {
            size_t bit = std::min_element(columnsSum.begin(), columnsSum.end(), [tableSize](size_t a, size_t b) {
                return std::abs(2 * (int64_t) a - (int64_t) tableSize) <
                       std::abs(2 * (int64_t) b - (int64_t) tableSize);
            }) - columnsSum.begin();
            return bit;
        }
    }

    size_t MaxDispersionBitSelector::operator()(const std::vector<std::filesystem::path> &fpTablePaths) const {
        std::vector<std::future<std::pair<std::vector<uint64_t>, uint64_t>>> tasks;
        tasks.resize(fpTablePaths.size());
        for (auto &filePath: fpTablePaths) {
            tasks.emplace_back(std::async(std::launch::async, findColumnsSumInFile, filePath));
        }
        std::vector<size_t> columnsSum(IndigoFingerprint::size(), 0);
        size_t tableSize = 0;
        for (auto &task: tasks) {
            const auto &[fileColumnsSum, fileTableSize] = task.get();
            mergeResults(columnsSum, tableSize, fileColumnsSum, fileTableSize);
        }
        return selectBitByColumnsSum(columnsSum, tableSize);
    }

    size_t MaxDispersionBitSelector::operator()(const std::filesystem::path &fpTablePath) const {
        auto [columnsSum, tableSize] = findColumnsSumInFile(fpTablePath);
        return selectBitByColumnsSum(columnsSum, tableSize);
    }
} // qtr