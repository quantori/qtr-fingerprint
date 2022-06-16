#pragma once

#include <gtest/gtest.h>

#include <filesystem>

namespace qtr
{
    class GTestEnvironment final : public testing::Environment
    {
    public:
        GTestEnvironment(const std::filesystem::path &bigDataDirPath);
        void SetUp() override;

    private:
        std::filesystem::path _bigDataDirPath;
    };
}