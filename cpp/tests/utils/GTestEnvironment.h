#pragma once

#include <gtest/gtest.h>

#include <filesystem>

namespace qtr
{
    class GTestEnvironment final : public testing::Environment
    {
    public:
        GTestEnvironment(std::filesystem::path dataDirPath, std::filesystem::path bigDataDirPath);
        void SetUp() override;

    private:
        std::filesystem::path _bigDataDirPath;
        std::filesystem::path _dataDirPath;
    };
}