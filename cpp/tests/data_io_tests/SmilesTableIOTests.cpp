#include <filesystem>

#include "gtest/gtest.h"

#include "smiles_table_io/SmilesTableReader.h"
#include "smiles_table_io/SmilesTableWriter.h"

#include "../utils/DataPathManager.h"


class SmilesTableIOTests : public ::testing::Test {
protected:
    void writeTmpTable(const std::vector<qtr::smiles_table_value_t> &values) {
        qtr::SmilesTableWriter writer(_tmpSmilesTableFilePath);
        writer << values;
    }

    std::vector<qtr::smiles_table_value_t> readTmpTable() {
        qtr::SmilesTableReader reader(_tmpSmilesTableFilePath);
        std::vector<qtr::smiles_table_value_t> bucket;
        reader >> bucket;
        return bucket;
    }

    void SetUp() override {
        _tmpSmilesTableFilePath = qtr::DataPathManager::getTmpDataDir() / "SmilesTableTmp";
    }

    void TearDown() override {
        EXPECT_EQ(expected, actual);
        std::filesystem::remove(_tmpSmilesTableFilePath);
    }

    std::filesystem::path _tmpSmilesTableFilePath;
    std::vector<qtr::smiles_table_value_t> expected, actual;
};


TEST_F(SmilesTableIOTests, EmptyFile) {
    writeTmpTable(expected);
    actual = readTmpTable();
}

TEST_F(SmilesTableIOTests, DefaultValue) {
    auto defaultValue = std::make_pair(0, std::string());
    expected = {defaultValue};
    writeTmpTable(expected);
    actual = readTmpTable();
}


TEST_F(SmilesTableIOTests, DefaultValues) {
    auto defaultValue = std::make_pair(0, std::string());
    expected = {defaultValue, defaultValue, defaultValue};
    writeTmpTable(expected);
    actual = readTmpTable();
}

TEST_F(SmilesTableIOTests, RandomValue) {
    std::string str1 = "C1=CC=C(C=C1)C=O";
    std::string str2 = "InChI=1S/C3H6O/c1-3(2)4/h1-2H3";
    std::string str3 = "C9H8O4";
    expected = {std::make_pair(0, str1), std::make_pair(1, str2), std::make_pair(2, str3)};
    writeTmpTable(expected);
    actual = readTmpTable();
}
