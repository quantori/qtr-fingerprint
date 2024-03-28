#include "ColumnsStatistic.h"

#include <future>
#include <random>

#include "data_io/fingerprint_table_io/FingerprintTableReader.h"

namespace qtr {

    namespace {
        ColumnsStatistic buildColumnsStatistic(const std::filesystem::path &filePath) {
            return ColumnsStatistic(filePath);
        }
    }

    ColumnsStatistic::ColumnsStatistic(size_t columnsCount) : _zerosCount(columnsCount, 0), _fingerprintsCount(0) {}

    ColumnsStatistic::ColumnsStatistic(const std::filesystem::path &filePath) : ColumnsStatistic() {
        collectStatistic(filePath);
    }

    ColumnsStatistic::ColumnsStatistic(const std::vector<std::filesystem::path> &filePaths) : ColumnsStatistic() {
        collectStatistic(filePaths);
    }

    size_t ColumnsStatistic::zeros(size_t i) const {
        return _zerosCount[i];
    }

    size_t ColumnsStatistic::ones(size_t i) const {
        return _fingerprintsCount - zeros(i);
    }

    size_t ColumnsStatistic::size() const {
        return _fingerprintsCount;
    }

    void ColumnsStatistic::collectStatistic(const std::filesystem::path &filePath) {
        FingerprintTableReader reader(filePath);
        for (const auto &[_, fingerprint]: reader) {
            if (_zerosCount.empty())
                _zerosCount.resize(fingerprint.size());
            else
                assert(_zerosCount.size() == fingerprint.size());
            for (size_t i = 0; i < columns(); i++) {
                _zerosCount[i] += fingerprint[i];
            }
            _fingerprintsCount++;
        }
    }

    void ColumnsStatistic::collectStatistic(const std::vector<std::filesystem::path> &filePaths) {
        std::vector<std::future<ColumnsStatistic>> tasks;
        for (auto &filePath: filePaths) {
            tasks.emplace_back(std::async(std::launch::async, buildColumnsStatistic, filePath));
        }
        for (auto &task: tasks) {
            operator+=(task.get());
        }
    }

    ColumnsStatistic &ColumnsStatistic::operator+=(const ColumnsStatistic &otherStatistic) {
        if (_zerosCount.empty()) {
            return *this = otherStatistic;
        }
        if (otherStatistic._zerosCount.empty()) {
            return *this;
        }
        assert(columns() == otherStatistic.columns());
        _fingerprintsCount += otherStatistic._fingerprintsCount;
        for (size_t i = 0; i < columns(); i++) {
            _zerosCount[i] += otherStatistic._zerosCount[i];
        }
        return *this;
    }

    size_t ColumnsStatistic::columns() const {
        return _zerosCount.size();
    }

    ColumnsStatistic::ColumnsStatistic() : ColumnsStatistic(0) {

    }


} // qtr