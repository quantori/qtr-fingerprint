#include "gtest/gtest.h"

#include "data_io/fingerprint_table_io/FingerprintTableWriter.h"
#include "data_io/fingerprint_table_io/FingerprintTableReader.h"

#include "../utils/DataPathManager.h"


class FingerprintTableIOTests : public ::testing::Test {
protected:
    void writeTmpTable(const std::vector<qtr::fingerprint_table_value_t> &values) {
        qtr::FingerprintTableWriter writer(_tmpFingerprintTableFilePath);
        writer << values;
    }

    std::vector<qtr::fingerprint_table_value_t> readTmpTable() {
        qtr::FingerprintTableReader reader(_tmpFingerprintTableFilePath);
        std::vector<qtr::fingerprint_table_value_t> table;
        reader >> table;
        return table;
    }

    void SetUp() override {
        _tmpFingerprintTableFilePath = qtr::DataPathManager::getTmpDataDir() / "FingerprintTableTmp";
    }

    void TearDown() override {
        EXPECT_EQ(expected, actual);
        std::filesystem::remove(_tmpFingerprintTableFilePath);
    }

    std::filesystem::path _tmpFingerprintTableFilePath;
    std::vector<qtr::fingerprint_table_value_t> expected, actual;
};

TEST_F(FingerprintTableIOTests, EmptyFile) {
    writeTmpTable(expected);
    actual = readTmpTable();
}

TEST_F(FingerprintTableIOTests, DefaultValue) {
    auto defaultValue = std::make_pair(0, qtr::IndigoFingerprint());
    expected = {defaultValue};
    writeTmpTable(expected);
    actual = readTmpTable();
}


TEST_F(FingerprintTableIOTests, DefaultValues) {
    auto defaultValue = std::make_pair(0, qtr::IndigoFingerprint());
    expected = {defaultValue, defaultValue, defaultValue};
    writeTmpTable(expected);
    actual = readTmpTable();
}

TEST_F(FingerprintTableIOTests, RandomValue) {
    auto fp1 = qtr::IndigoFingerprint(), fp2 = qtr::IndigoFingerprint(), fp3 = qtr::IndigoFingerprint();
    fp1[0] = fp1[1] = fp1[10] = fp1[20] = fp1[qtr::IndigoFingerprint::size() - 1] = true;
    fp2[1] = fp2[3] = fp2[4] = true;
    fp3[100] = fp3[200] = fp3[400] = true;
    expected = {std::make_pair(0, fp1), std::make_pair(1, fp2), std::make_pair(2, fp3)};
    writeTmpTable(expected);
    actual = readTmpTable();
}