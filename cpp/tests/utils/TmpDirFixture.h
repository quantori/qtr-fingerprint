#pragma once

#include "gtest/gtest.h"

#include "DataPathManager.h"


class TmpDirFixture : public ::testing::Test {
public:
    void SetUp() override;

    void TearDown() override;

    const std::filesystem::path& getTmpDir() const;
private:
    std::filesystem::path _dir;

};


