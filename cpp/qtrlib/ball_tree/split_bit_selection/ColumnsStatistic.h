#pragma once

#include <vector>
#include <filesystem>

namespace qtr {

    class ColumnsStatistic {
    public:
        ColumnsStatistic();

        explicit ColumnsStatistic(const std::filesystem::path &filePath);

        explicit ColumnsStatistic(const std::vector<std::filesystem::path> &filePaths);

        [[nodiscard]] size_t zeros(size_t i) const;

        [[nodiscard]] size_t ones(size_t i) const;

        [[nodiscard]] size_t size() const;

        void collectStatistic(const std::filesystem::path &filePath);

        void collectStatistic(const std::vector<std::filesystem::path> &filePaths);

        ColumnsStatistic &operator+=(const ColumnsStatistic &otherStatistic);

    private:
        std::vector<size_t> _zerosCount;
        size_t _fingerprintsCount;
    };

} // qtr
