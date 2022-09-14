#include <filesystem>

#include "gtest/gtest.h"
#include "../utils/DataPathManager.h"

#include "ColumnsIO.h"


class ColumnsIOTests : public ::testing::Test {
protected:
    void writeTmpColumns(const std::vector<size_t> &values) {
        qtr::ColumnsWriter writer(_tmpColumnsFilePath);
        writer.write(values);
    }

    std::vector<size_t> readTmpColumns() {
        qtr::ColumnsReader reader(_tmpColumnsFilePath);
        return reader.readAll();
    }

    void SetUp() override {
        _tmpColumnsFilePath = qtr::DataPathManager::getTmpDataDir() / "ColumnsTmp";
    }

    void TearDown() override {
        EXPECT_EQ(expected, actual);
        std::filesystem::remove(_tmpColumnsFilePath);
    }

    std::filesystem::path _tmpColumnsFilePath;
    std::vector<size_t> expected, actual;
};



TEST_F(ColumnsIOTests, EmptyFile) {
    writeTmpColumns(expected);
    actual = readTmpColumns();
}

TEST_F(ColumnsIOTests, Value) {
    expected = {1};
    writeTmpColumns(expected);
    actual = readTmpColumns();
}

TEST_F(ColumnsIOTests, Values) {
    expected = {5, 6, 4, 1, 3, 0, 2};
    writeTmpColumns(expected);
    actual = readTmpColumns();
}
