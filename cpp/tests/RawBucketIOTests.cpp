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

TEST(RAwBucketIOTest, RealValue) {
    // TODO test with real value
}

TEST(RAwBucketIOTest, RealValues) {
    // TODO test with real values
}
