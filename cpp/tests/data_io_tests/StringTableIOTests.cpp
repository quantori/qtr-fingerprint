#include <filesystem>

#include "gtest/gtest.h"

#include "string_table_io/StringTableReader.h"
#include "string_table_io/StringTableWriter.h"

#include "../utils/DataPathManager.h"


class StringTableIOTests : public ::testing::Test {
protected:
    void writeTmpTable(const std::vector<qtr::string_table_value_t> &values) {
        qtr::StringTableWriter writer(_tmpTableFilePath);
        writer << values;
    }

    std::vector<qtr::string_table_value_t> readTmpTable() {
        qtr::StringTableReader reader(_tmpTableFilePath);
        std::vector<qtr::string_table_value_t> bucket;
        reader >> bucket;
        return bucket;
    }

    void SetUp() override {
        _tmpTableFilePath = qtr::DataPathManager::getTmpDataDir() / "StringTableTmp";
    }

    void TearDown() override {
        EXPECT_EQ(expected, actual);
        std::filesystem::remove(_tmpTableFilePath);
    }

    std::filesystem::path _tmpTableFilePath;
    std::vector<qtr::string_table_value_t> expected, actual;
};


TEST_F(StringTableIOTests, EmptyFile) {
    writeTmpTable(expected);
    actual = readTmpTable();
}

TEST_F(StringTableIOTests, DefaultValue) {
    auto defaultValue = std::make_pair(0, std::string());
    expected = {defaultValue};
    writeTmpTable(expected);
    actual = readTmpTable();
}


TEST_F(StringTableIOTests, DefaultValues) {
    auto defaultValue = std::make_pair(0, std::string());
    expected = {defaultValue, defaultValue, defaultValue};
    writeTmpTable(expected);
    actual = readTmpTable();
}

TEST_F(StringTableIOTests, RandomValue) {
    std::string str1 = "C1=CC=C(C=C1)C=O";
    std::string str2 = "InChI=1S/C3H6O/c1-3(2)4/h1-2H3";
    std::string str3 = "C9H8O4";
    expected = {std::make_pair(0, str1), std::make_pair(1, str2), std::make_pair(2, str3)};
    writeTmpTable(expected);
    actual = readTmpTable();
}
