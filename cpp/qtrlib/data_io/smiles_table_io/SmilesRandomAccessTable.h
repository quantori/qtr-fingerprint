#pragma once

#include <vector>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>

namespace qtr {

    class SmilesRandomAccessTable {
    public:
        static inline std::filesystem::path getSmilesStoragePath(const std::filesystem::path &tableDir) {
            return tableDir / "storage";
        }

        static inline std::filesystem::path getSeeksPath(const std::filesystem::path &tableDir) {
            return tableDir / "seeks.txt";
        }

        SmilesRandomAccessTable(const std::filesystem::path& smilesTablePath, const std::filesystem::path &tableDir);

        explicit SmilesRandomAccessTable(const std::filesystem::path &tableDir);

        std::string operator[](size_t i);

        void toVector(std::vector<std::string>& vector);

    private:
        std::vector<uint64_t> _seeks;
        std::ifstream _smilesFile;
    };

}