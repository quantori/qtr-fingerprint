#pragma once

#include <vector>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>

namespace qtr {

    class StringRandomAccessTable {
    public:
        static inline std::filesystem::path getStoragePath(const std::filesystem::path &tableDir) {
            return tableDir / "storage";
        }

        static inline std::filesystem::path getSeeksPath(const std::filesystem::path &tableDir) {
            return tableDir / "seeks.txt";
        }

        StringRandomAccessTable(const std::filesystem::path &tablePath, const std::filesystem::path &tableDir);

        explicit StringRandomAccessTable(const std::filesystem::path &tableDir);

        std::string operator[](size_t i);

    private:
        std::vector<uint64_t> _seeks;
        std::ifstream _tableFile;
    };

}