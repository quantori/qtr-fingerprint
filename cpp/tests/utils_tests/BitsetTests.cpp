#include "gtest/gtest.h"

#include "../utils/DataPathManager.h"

#include "Fingerprint.h"
#include "io/BufferedReader.h"
#include "io/BufferedWriter.h"

class BitsetTests : public ::testing::Test {
public:
    static const size_t BitsetSize = 17;
    using BitsetStorageType = unsigned char;
    qtr::Bitset<BitsetSize, BitsetStorageType> bitset;

    void changeValAndCheck(size_t pos, bool newVal) {
        std::vector<bool> oldValues(BitsetSize);
        for (size_t i = 0; i < BitsetSize; i++) {
            oldValues[i] = bitset[i];
        }
        bitset[pos] = newVal;
        for (size_t i = 0; i < BitsetSize; i++) {
            if (i == pos) {
                EXPECT_EQ(bitset[i], newVal);
            } else {
                EXPECT_EQ(bitset[i], oldValues[i]);
            }
        }
    }

};

TEST_F(BitsetTests, DefaultInitializationTest) {
    size_t bitsetSize = BitsetSize; // copy size to use it on #define EXPECT_EQ
    EXPECT_EQ(bitsetSize, bitset.size());
    for (size_t i = 0; i < BitsetSize; i++) {
        EXPECT_FALSE(bitset[i]);
    }
}


TEST_F(BitsetTests, OperatorSquareBracketsTests) {
    changeValAndCheck(0, true);
    changeValAndCheck(0, true);
    changeValAndCheck(0, false);
    changeValAndCheck(0, false);
    for (size_t i = 0; i < BitsetSize; i++) {
        changeValAndCheck(i, true);
    }
    for (size_t i = 0; i < BitsetSize; i++) {
        changeValAndCheck(BitsetSize - i, false);
    }
}

TEST_F(BitsetTests, CopyAsignmentTest) {
    for (size_t i = 0; i < 3; i++) {
        bitset[i] = true;
    }
    auto bitset2 = bitset;
    for (size_t i = 3; i < 5; i++) {
        bitset[i] = true;
    }
    for (size_t i = 0; i < BitsetSize; i++) {
        EXPECT_EQ(bitset2[i], i < 3);
    }
}

TEST_F(BitsetTests, OperatorEqualEqualTests) {
    EXPECT_TRUE(bitset == bitset);
    for (size_t i = 0; i < 3; i++) {
        bitset[i] = true;
    }
    auto bitset2 = bitset;
    EXPECT_TRUE(bitset2 == bitset);
    EXPECT_TRUE(bitset2 == bitset2);
    bitset[4] = true;
    EXPECT_FALSE(bitset2 == bitset);
    bitset2[4] = true;
    EXPECT_TRUE(bitset == bitset2);
}

TEST_F(BitsetTests, OperatorOrEqualTests) {
    auto bitset2 = bitset;
    for (size_t i = 6; i < 10; i++)
        bitset2[i] = true;
    for (size_t i = 2; i < 4; i++)
        bitset[i] = true;
    bitset |= bitset2;
    for (size_t i = 0; i < BitsetSize; i++) {
        EXPECT_EQ(bitset[i], (2 <= i && i < 4) || (6 <= i && i < 10));
    }
    bitset2 |= bitset2;
    for (size_t i = 0; i < BitsetSize; i++) {
        EXPECT_EQ(bitset2[i], 6 <= i && i < 10);
    }
    bitset2 |= bitset;
    for (size_t i = 0; i < BitsetSize; i++) {
        EXPECT_EQ(bitset2[i], (2 <= i && i < 4) || (6 <= i && i < 10));
    }
}

TEST_F(BitsetTests, OperatorOrTests) {
    auto bitset2 = bitset;
    for (size_t i = 6; i < 10; i++)
        bitset2[i] = true;
    for (size_t i = 2; i < 4; i++)
        bitset[i] = true;
    for (size_t i = 0; i < BitsetSize; i++) {
        EXPECT_EQ((bitset | bitset2)[i], (2 <= i && i < 4) || (6 <= i && i < 10));
    }
    for (size_t i = 0; i < BitsetSize; i++) {
        EXPECT_EQ((bitset2 | bitset)[i], (2 <= i && i < 4) || (6 <= i && i < 10));
    }
    for (size_t i = 0; i < BitsetSize; i++) {
        EXPECT_EQ((bitset2 | bitset2)[i], 6 <= i && i < 10);
    }
}

TEST_F(BitsetTests, OperatorLowerEqualTest) {
    EXPECT_TRUE(bitset <= bitset);
    for (size_t i = 0; i < 3; i++) {
        bitset[i] = true;
    }
    auto bitset2 = bitset;
    EXPECT_TRUE(bitset2 <= bitset);
    bitset[8] = true;
    EXPECT_TRUE(bitset2 <= bitset);
    EXPECT_FALSE(bitset <= bitset2);
    bitset2[7] = true;
    EXPECT_FALSE(bitset2 <= bitset);
    EXPECT_FALSE(bitset <= bitset2);
    bitset2[8] = true;
    EXPECT_FALSE(bitset2 <= bitset);
    EXPECT_TRUE(bitset <= bitset2);
    bitset[7] = true;
    EXPECT_TRUE(bitset2 <= bitset);
    EXPECT_TRUE(bitset <= bitset2);
}

TEST_F(BitsetTests, StdIODumpLoadTest) {
    auto filePath = qtr::DataPathManager::getTmpDataDir() / "BitsetTmp";
    for (size_t i = 6; i < 10; i++)
        bitset[i] = true;
    for (size_t i = 2; i < 4; i++)
        bitset[i] = true;
    auto bitset2 = bitset;
    {
        std::ofstream writer(filePath);
        bitset.dump(writer);
    }
    {
        std::ifstream reader(filePath);
        bitset.reset();
        bitset.load(reader);
        EXPECT_EQ(bitset2, bitset);
    }
    std::filesystem::remove(filePath);
}

TEST_F(BitsetTests, QtrIODumpLoadTest) {
    auto filePath = qtr::DataPathManager::getTmpDataDir() / "BitsetTmp";
    for (size_t i = 6; i < 10; i++)
        bitset[i] = true;
    for (size_t i = 2; i < 4; i++)
        bitset[i] = true;
    auto bitset2 = bitset;
    {
        qtr::BufferedWriter<1024> writer(filePath);
        bitset.dump(writer);
    }
    {
        qtr::BufferedReader<1024> reader(filePath);
        bitset.reset();
        bitset.load(reader);
        EXPECT_EQ(bitset2, bitset);
    }
    std::filesystem::remove(filePath);
}

