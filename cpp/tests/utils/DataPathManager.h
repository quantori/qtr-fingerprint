#pragma once

#include <filesystem>

namespace qtr
{
    class DataPathManager
    {
    public:
        static std::filesystem::path getDataDir();
        static std::filesystem::path getBigDataDir();
        static void init(const std::filesystem::path &bigDataDirPath);

    private:
        inline static std::filesystem::path _bigDataDirPath;
    };
}