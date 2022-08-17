#pragma once

#include <filesystem>

namespace qtr
{
    class DataPathManager
    {
    public:
        static std::filesystem::path getDataDir();
        static std::filesystem::path getBigDataDir();
        static std::filesystem::path getTmpDataDir();
        static void init(const std::filesystem::path &dataDirPath, const std::filesystem::path &bigDataDirPath, const std::filesystem::path& tmpDataDirPath);

    private:
        inline static std::filesystem::path _dataDirPath;
        inline static std::filesystem::path _bigDataDirPath;
        inline static std::filesystem::path _tmpDataDirPath;
    };
}
