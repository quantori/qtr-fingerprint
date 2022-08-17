#pragma once

#include <gtest/gtest.h>

#include <filesystem>

namespace qtr {
    class GTestEnvironment final : public testing::Environment {
    public:
        GTestEnvironment(std::filesystem::path dataDirPath, std::filesystem::path bigDataDirPath,
                         std::filesystem::path tmpDataDirPath);

        void SetUp() override;

    private:
        std::filesystem::path _dataDirPath;
        std::filesystem::path _bigDataDirPath;
        std::filesystem::path _tmpDataDirPath;
    };
}