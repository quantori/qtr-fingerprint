#include <filesystem>

#include "gtest/gtest.h"

#include "raw_bucket_io/RawBucketReader.h"
#include "raw_bucket_io/RawBucketWriter.h"

#include "../utils/DataPathManager.h"


class RawBucketIOTests : public ::testing::Test {
protected:
    void writeTmpBucket(const std::vector<qtr::raw_bucket_value_t> &values) {
        qtr::RawBucketWriter writer(_tmpRawBucketFilePath);
        writer.write(values);
    }

    std::vector<qtr::raw_bucket_value_t> readTmpBucket() {
        qtr::RawBucketReader reader(_tmpRawBucketFilePath);
        return reader.readAll();
    }

    void SetUp() override {
        _tmpRawBucketFilePath = qtr::DataPathManager::getTmpDataDir() / "RawBucketTmp";
    }

    void TearDown() override {
        EXPECT_EQ(expected, actual);
        std::filesystem::remove(_tmpRawBucketFilePath);
    }

    std::filesystem::path _tmpRawBucketFilePath;
    std::vector<qtr::raw_bucket_value_t> expected, actual;
};


TEST_F(RawBucketIOTests, EmptyFile) {
    writeTmpBucket(expected);
    actual = readTmpBucket();
}

TEST_F(RawBucketIOTests, DefaultValue) {
    auto defaultValue = std::make_pair(std::string(), qtr::IndigoFingerprint());
    expected = {defaultValue};
    writeTmpBucket(expected);
    actual = readTmpBucket();
}


TEST_F(RawBucketIOTests, DefaultValues) {
    auto defaultValue = std::make_pair(std::string(), qtr::IndigoFingerprint());
    expected = {defaultValue, defaultValue, defaultValue};
    writeTmpBucket(expected);
    actual = readTmpBucket();
}

TEST_F(RawBucketIOTests, RandomValue) {
    auto fp1 = qtr::IndigoFingerprint(), fp2 = qtr::IndigoFingerprint(), fp3 = qtr::IndigoFingerprint();
    fp1[0] = fp1[1] = fp1[10] = fp1[20] = fp1[qtr::IndigoFingerprint::size() - 1] = true;
    fp2[1] = fp2[3] = fp2[4] = true;
    fp3[100] = fp3[200] = fp3[400] = true;
    std::string str1 = "C1=CC=C(C=C1)C=O";
    std::string str2 = "InChI=1S/C3H6O/c1-3(2)4/h1-2H3";
    std::string str3 = "C9H8O4";
    expected = {std::make_pair(str1, fp1), std::make_pair(str2, fp2), std::make_pair(str3, fp3)};
    writeTmpBucket(expected);
    actual = readTmpBucket();
}
