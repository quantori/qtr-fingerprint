#include "gtest/gtest.h"

#include "SplitterTree.h"
#include "SplitterTreeUtils.h"
#include "FingerprintTable.h"
#include "utils/DataPathManager.h"
#include "raw_bucket_io/RawBucketReader.h"
#include "raw_bucket_io/RawBucketWriter.h"

using namespace qtr;

class SplitterTreeTests : public ::testing::Test {
protected:
    void SetUp() override {
        fingerprints.resize(5);
        // fingerprints:
        // 0: 1000...
        // 1: 1010...
        // 2: 1100...
        // 3: 1110...
        // 4: 1100...
        fingerprints[0][0] = true;
        fingerprints[1][0] = fingerprints[1][2] = true;
        fingerprints[2][0] = fingerprints[2][1] = true;
        fingerprints[3][0] = fingerprints[3][1] = fingerprints[3][2] = true;
        fingerprints[4][0] = fingerprints[4][1] = true;

        rawBucketValues.resize(5); // SMILES have been chosen randomly
        rawBucketValues[0] = {"C1=C(C=C(C(=C1O)O)O)C2=[O+]C3=CC(=CC(=C3C=C2O)O)O.[Cl-]", fingerprints[0]};
        rawBucketValues[1] = {"C1=CC=C(C=C1)C=O", fingerprints[1]};
        rawBucketValues[2] = {"C=C", fingerprints[2]};
        rawBucketValues[3] = {"[C-]#[O+]", fingerprints[3]};
        rawBucketValues[4] = {"CC1CCC2=C3N1C=C(C(=O)C3=CC(=C2)F)C(=O)O", fingerprints[4]};


        // splitter tree's buckets in nodes
        //        0
        //    1        2
        //  3   4    5   6

        rawBuckets.resize(7);
        rawBuckets[0] = rawBucketValues; // root bucket, split bit = 1, left child = 1, right child = 2
        rawBuckets[1] = {rawBucketValues[0], rawBucketValues[1]}; // split bit = 2, left child = 3, right child = 4
        rawBuckets[2] = {rawBucketValues[2], rawBucketValues[3],
                         rawBucketValues[4]}; // split bit = 2, left child = 5, right child = 6
        rawBuckets[3] = {rawBucketValues[0]}; // leaf node
        rawBuckets[4] = {rawBucketValues[1]}; // leaf node
        rawBuckets[5] = {rawBucketValues[2], rawBucketValues[4]}; // leaf node
        rawBuckets[6] = {rawBucketValues[3]}; // leaf node

        splitBits = {1, 2, 2, (uint64_t) -1, (uint64_t) -1, (uint64_t) -1, (uint64_t) -1};

        children.resize(7);
        children[0] = {1, 2};
        children[1] = {3, 4};
        children[2] = {5, 6};
        children[3] = {(uint64_t) -1, (uint64_t) -1};
        children[4] = {(uint64_t) -1, (uint64_t) -1};
        children[5] = {(uint64_t) -1, (uint64_t) -1};
        children[6] = {(uint64_t) -1, (uint64_t) -1};

        std::filesystem::create_directory(rawBucketsDirPath);
    }

    std::filesystem::path rawBucketPath(size_t i) {
        std::filesystem::create_directory(rawBucketsDirPath / std::to_string(i));
        return rawBucketsDirPath / std::to_string(i) / "data.rb";
    }

    static void dumpBucket(const std::filesystem::path &bucketPath, const std::vector<raw_bucket_value_t> &bucket) {
        RawBucketWriter writer(bucketPath);
        writer << bucket;
    }

    void dumpAllBuckets() {
        for (size_t i = 0; i < rawBuckets.size(); i++) {
            dumpBucket(rawBucketPath(i), rawBuckets[i]);
        }
    }

    static std::vector<raw_bucket_value_t> loadBucket(const std::filesystem::path &bucketPath) {
        RawBucketReader reader(bucketPath);
        std::vector<raw_bucket_value_t> bucket;
        reader >> bucket;
        return bucket;
    }

    static void
    checkBucketsEqual(const std::vector<raw_bucket_value_t> &bucket1, const std::vector<raw_bucket_value_t> &bucket2) {
        EXPECT_TRUE(std::is_permutation(bucket1.begin(), bucket1.end(), bucket2.begin(), bucket2.end()));
    }

    void TearDown() override {
        std::filesystem::remove_all(rawBucketsDirPath);
    }

    IndigoFingerprintTable fingerprints;
    std::vector<raw_bucket_value_t> rawBucketValues;
    std::vector<std::vector<raw_bucket_value_t>> rawBuckets;
    std::vector<uint64_t> splitBits;
    std::vector<std::pair<uint64_t, uint64_t>> children;
    std::filesystem::path rawBucketsDirPath = DataPathManager::getTmpDataDir() / "rawBucketsTmp";
};


TEST_F(SplitterTreeTests, FindBestBitToSplitTest) {
    dumpAllBuckets();
    for (size_t i = 0; i < 3; i++) {
        uint64_t splitBit = findBestBitToSplit({rawBucketPath(i)}, false);
        EXPECT_EQ(splitBit, splitBits[i]);
    }
}

TEST_F(SplitterTreeTests, SplitRawBucketByBitTest) {
    dumpAllBuckets();
    for (size_t i = 0; i < 3; i++) {
        auto zerosBucketPath = rawBucketsDirPath / "zerosBucketTmp";
        auto onesBucketPath = rawBucketsDirPath / "onesBucketTmp";
        splitRawBucketByBitNotParallel({rawBucketPath(i)}, splitBits[i], zerosBucketPath, onesBucketPath);
        auto actualZerosBucket = loadBucket(zerosBucketPath);
        auto actualOnesBucket = loadBucket(onesBucketPath);
        auto expectedZerosBucket = rawBuckets[i * 2 + 1];
        auto expectedOnesBucket = rawBuckets[i * 2 + 2];
        checkBucketsEqual(actualZerosBucket, expectedZerosBucket);
        checkBucketsEqual(actualOnesBucket, expectedOnesBucket);
    }
}


TEST_F(SplitterTreeTests, BuildSmallDebthNotParallelTest) {
    dumpBucket(rawBucketPath(0), rawBuckets[0]);
    SplitterTree tree({rawBucketsDirPath});
    tree.build(1, 2, 0);
    EXPECT_EQ(tree.size(), 3);
    for (size_t i = 1; i <= 2; i++) {
        auto actualBucket = loadBucket(rawBucketPath(i));
        auto &expectedBucket = rawBuckets[i];
        checkBucketsEqual(actualBucket, expectedBucket);
    }
}

TEST_F(SplitterTreeTests, BuildSmallDebthParallelTest) {
    dumpBucket(rawBucketPath(0), rawBuckets[0]);
    SplitterTree tree({rawBucketsDirPath});
    tree.build(1, 1, 0);
    EXPECT_EQ(tree.size(), 3);
    for (size_t i = 1; i <= 2; i++) {
        auto actualBucket = loadBucket(rawBucketPath(i));
        auto &expectedBucket = rawBuckets[i];
        checkBucketsEqual(actualBucket, expectedBucket);
    }
}

TEST_F(SplitterTreeTests, BuildNotParallelTest) {
    dumpBucket(rawBucketPath(0), rawBuckets[0]);
    SplitterTree tree({rawBucketsDirPath});
    tree.build(2, 1, 0);
    EXPECT_EQ(tree.size(), 7);
    std::ofstream out(rawBucketsDirPath);
    tree.dump(out);
    for (size_t i = 3; i <= 6; i++) {
        auto actualBucket = loadBucket(rawBucketPath(i));
        auto &expectedBucket = rawBuckets[i];
        checkBucketsEqual(actualBucket, expectedBucket);
    }
}

TEST_F(SplitterTreeTests, DumpTest) {
    dumpBucket(rawBucketPath(0), rawBuckets[0]);
    SplitterTree tree({rawBucketsDirPath});
    tree.build(2, 1, 1);
    auto filePath = DataPathManager::getTmpDataDir() / "splitterTreeTmp";
    {
        std::ofstream out(filePath);
        tree.dump(out);
    }
    {
        std::ifstream in(filePath);

        auto getNum = [&in]() {
            uint64_t num;
            in.read((char *) &num, sizeof num);
            return num;
        };

        uint64_t treeSize = getNum();
        EXPECT_EQ(treeSize, 7);
        for (size_t i = 0; i < 7; i++) {
            uint64_t nodeId = getNum();
            uint64_t splitBit = getNum();
            uint64_t leftChild = getNum();
            uint64_t rightChild = getNum();
            EXPECT_LT(nodeId, 7);
            EXPECT_EQ(splitBit, splitBits[nodeId]);
            EXPECT_EQ(leftChild, children[nodeId].first);
            EXPECT_EQ(rightChild, children[nodeId].second);
        }
    }
    std::filesystem::remove(filePath);
}
