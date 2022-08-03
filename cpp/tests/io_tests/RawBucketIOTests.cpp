#include <filesystem>

#include "gtest/gtest.h"

#include "RawBucketsIO.h"

//static const std::filesystem::path tmpBucketPath = std::filesystem::temp_directory_path() / "rawBucketTmp";
static const std::filesystem::path tmpBucketPath = "/home/Vsevolod.Vaskin/qtr-fingerprint/data/tests/rawBucketTmp";

static void writeTmpBucket(const std::vector<qtr::raw_bucket_value_t> &values) {
    qtr::RawBucketWriter writer(tmpBucketPath);
    writer.write(values);
}

static std::vector<qtr::raw_bucket_value_t> readTmpBucket() {
    qtr::RawBucketReader reader(tmpBucketPath);
    return reader.readAll();
}

TEST(RawBucketIOTest, EmptyFile) {
    std::vector<qtr::raw_bucket_value_t> expected;
    writeTmpBucket(expected);
    auto actual = readTmpBucket();
    EXPECT_EQ(expected, actual);
}

TEST(RAwBucketIOTest, DefaultValue) {
    auto defaultValue = std::make_pair(std::string(), qtr::IndigoFingerprint());
    std::vector<qtr::raw_bucket_value_t> expected = {defaultValue};
    writeTmpBucket(expected);
    auto actual = readTmpBucket();
    EXPECT_EQ(expected, actual);
}


TEST(RAwBucketIOTest, DefaultValues) {
    auto defaultValue = std::make_pair(std::string(), qtr::IndigoFingerprint());
    std::vector<qtr::raw_bucket_value_t> expected = {defaultValue, defaultValue, defaultValue};
    writeTmpBucket(expected);
    auto actual = readTmpBucket();
    EXPECT_EQ(expected, actual);
}

TEST(RAwBucketIOTest, RandomValue) {
    auto fp1 = qtr::IndigoFingerprint(), fp2 = qtr::IndigoFingerprint(), fp3 = qtr::IndigoFingerprint();
    fp1[0] = fp1[1] = fp1[10] = fp1[20] = fp1[fp1.size() - 1] = true;
    fp2[1] = fp2[3] = fp2[4] = true;
    fp3[100] = fp3[200] = fp3[400] = true;
    std::string str1 = "C1=CC=C(C=C1)C=O";
    std::string str2 = "InChI=1S/C3H6O/c1-3(2)4/h1-2H3";
    std::string str3 = "C9H8O4";
    std::vector<qtr::raw_bucket_value_t> expected = {std::make_pair(str1, fp1),
                                                     std::make_pair(str2, fp2),
                                                     std::make_pair(str3, fp3)};
    writeTmpBucket(expected);
    auto actual = readTmpBucket();
    EXPECT_EQ(expected, actual);
}
