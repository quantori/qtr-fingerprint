#include <filesystem>

#include "gtest/gtest.h"
#include "../utils/DataPathManager.h"

#include "indexes_io/IndexesReader.h"
#include "indexes_io/IndexesWriter.h"


class IndexesIOTests : public ::testing::Test {
protected:
    void writeTmpIndexes(const std::vector<uint64_t> &values) {
        qtr::IndexesWriter writer(_tmpIndexesFilePath);
        writer << values;
    }

    std::vector<uint64_t> readTmpIndexes() {
        qtr::IndexesReader reader(_tmpIndexesFilePath);
        std::vector<qtr::IndexesReader::ReadValue> result;
        reader >> result;
        return result;
    }

    void SetUp() override {
        _tmpIndexesFilePath = qtr::DataPathManager::getTmpDataDir() / "IndexesTmp";
    }

    void TearDown() override {
        EXPECT_EQ(expected, actual);
        std::filesystem::remove(_tmpIndexesFilePath);
    }

    std::filesystem::path _tmpIndexesFilePath;
    std::vector<uint64_t> expected, actual;
};



TEST_F(IndexesIOTests, EmptyFile) {
    writeTmpIndexes(expected);
    actual = readTmpIndexes();
}

TEST_F(IndexesIOTests, Value) {
    expected = {1};
    writeTmpIndexes(expected);
    actual = readTmpIndexes();
}

TEST_F(IndexesIOTests, Values) {
    expected = {5, 6, 4, 1, 3, 0, 2};
    writeTmpIndexes(expected);
    actual = readTmpIndexes();
}
