#include <gtest/gtest.h>

#include "Utils.h"
#include "../utils/TmpDirFixture.h"

#include <bitset>
#include <fstream>
#include <vector>
#include <set>


static const int VALUE_SIZE_BITS = qtr::fromBytesToBits(sizeof(uint64_t));

TEST(getBit, ALL_ZEROES) {
    std::bitset<VALUE_SIZE_BITS> curr = 0;
    for (int i = 0; i < VALUE_SIZE_BITS; ++i) {
        EXPECT_EQ(curr[i], 0);
        EXPECT_EQ(curr[i], qtr::getBit(curr.to_ullong(), i));
    }
}

TEST(getBit, ALL_ONES) {
    std::bitset<VALUE_SIZE_BITS> curr = -1;
    for (int i = 0; i < VALUE_SIZE_BITS; ++i) {
        EXPECT_EQ(curr[i], 1);
        EXPECT_EQ(curr[i], qtr::getBit(curr.to_ullong(), i));
    }
}

TEST(getBit, RANDOM) {
    std::bitset<VALUE_SIZE_BITS> curr = 17432891238;
    for (int i = 0; i < VALUE_SIZE_BITS; ++i) {
        EXPECT_EQ(curr[i], qtr::getBit(curr.to_ullong(), i));
    }
}

TEST(fromBytesToBits, COMMON) {
    EXPECT_EQ(qtr::fromBytesToBits(0), 0);
    EXPECT_EQ(qtr::fromBytesToBits(1), 8);
    EXPECT_EQ(qtr::fromBytesToBits(1000), 8000);
}

TEST(divideIntegersCeil, COMMON) {
    EXPECT_EQ(qtr::divideIntegersCeil(0, 5), 0);
    EXPECT_EQ(qtr::divideIntegersCeil(1, 2), 1);
    EXPECT_EQ(qtr::divideIntegersCeil(1, 3), 1);
    EXPECT_EQ(qtr::divideIntegersCeil(3, 3), 1);
    EXPECT_EQ(qtr::divideIntegersCeil(4, 3), 2);
    EXPECT_EQ(qtr::divideIntegersCeil(15, 4), 4);
}

TEST_F(TmpDirFixture, findFilesTests) {
    std::filesystem::path tmpDir = getTmpDir();

    std::filesystem::path dirPath = tmpDir / "dir_name";
    std::filesystem::path fileNoExtensionPath = tmpDir / "filename";
    std::filesystem::path fileTxt1Path = tmpDir / "file1.txt";
    std::filesystem::path fileTxt2Path = tmpDir / "file2.txt";
    std::filesystem::path filePngPath = tmpDir / "file.png";
    std::filesystem::path fileEmptyExtensionPath = tmpDir / "file.";

    {
        std::filesystem::create_directory(dirPath);
        std::ofstream fileNoExtension(fileNoExtensionPath);
        std::ofstream fileTxt(fileTxt1Path);
        std::ofstream file2Txt(fileTxt2Path);
        std::ofstream filePng(filePngPath);
        std::ofstream fileEmptyExtension(fileEmptyExtensionPath);
    }

    EXPECT_EQ(0, qtr::findFiles(tmpDir, ".jpg").size());
    EXPECT_EQ(std::vector{tmpDir / "file.png"}, qtr::findFiles(tmpDir, ".png"));
    EXPECT_EQ(std::vector{tmpDir / "file.png"}, qtr::findFiles(tmpDir, "png"));

    std::vector<std::filesystem::path> txtFiles = {fileTxt1Path, fileTxt2Path};
    std::vector<std::filesystem::path> foundTxt1 = qtr::findFiles(tmpDir, ".txt");
    std::vector<std::filesystem::path> foundTxt2 = qtr::findFiles(tmpDir, "txt");
    EXPECT_TRUE(std::is_permutation(txtFiles.begin(), txtFiles.end(), foundTxt1.begin(), foundTxt1.end()));
    EXPECT_TRUE(std::is_permutation(txtFiles.begin(), txtFiles.end(), foundTxt2.begin(), foundTxt2.end()));

    std::set<std::filesystem::path> noExtension = {dirPath, fileNoExtensionPath};
    std::vector<std::filesystem::path> foundNoExtensionVector = qtr::findFiles(tmpDir, "");
    std::set<std::filesystem::path> foundNoExtension(foundNoExtensionVector.begin(), foundNoExtensionVector.end());
    EXPECT_EQ(noExtension, foundNoExtension);

    EXPECT_EQ(std::vector{fileEmptyExtensionPath}, qtr::findFiles(tmpDir, "."));
}

TEST(endsWithTest, COMMON) {
    EXPECT_TRUE(qtr::endsWith("abcd", "cd"));
    EXPECT_TRUE(qtr::endsWith("abcd", "abcd"));
    EXPECT_TRUE(qtr::endsWith("", ""));
    EXPECT_TRUE(qtr::endsWith("abc", ""));
    EXPECT_FALSE(qtr::endsWith("abcd", "ab"));
    EXPECT_FALSE(qtr::endsWith("", "a"));
    EXPECT_FALSE(qtr::endsWith("ab", "c"));
}

TEST(chexToIntTest, COMMON) {
    EXPECT_EQ(qtr::chexToInt('0'), 0);
    EXPECT_EQ(qtr::chexToInt('1'), 1);
    EXPECT_EQ(qtr::chexToInt('2'), 2);
    EXPECT_EQ(qtr::chexToInt('3'), 3);
    EXPECT_EQ(qtr::chexToInt('4'), 4);
    EXPECT_EQ(qtr::chexToInt('5'), 5);
    EXPECT_EQ(qtr::chexToInt('6'), 6);
    EXPECT_EQ(qtr::chexToInt('7'), 7);
    EXPECT_EQ(qtr::chexToInt('8'), 8);
    EXPECT_EQ(qtr::chexToInt('9'), 9);
    EXPECT_EQ(qtr::chexToInt('a'), 10);
    EXPECT_EQ(qtr::chexToInt('b'), 11);
    EXPECT_EQ(qtr::chexToInt('c'), 12);
    EXPECT_EQ(qtr::chexToInt('d'), 13);
    EXPECT_EQ(qtr::chexToInt('e'), 14);
    EXPECT_EQ(qtr::chexToInt('f'), 15);
}

TEST(lowerOrderBitsTest, COMMON) {
    EXPECT_EQ(qtr::lowerOrderBits(0, 1), 0);
    EXPECT_EQ(qtr::lowerOrderBits(123, 0), 0);
    EXPECT_EQ(qtr::lowerOrderBits(0b1010, 2), 0b10);
    EXPECT_EQ(qtr::lowerOrderBits(0b10111, 4), 0b111);
    EXPECT_EQ(qtr::lowerOrderBits(0b10111, 10), 0b10111);
    EXPECT_EQ(qtr::lowerOrderBits(0xffffffffull, 32), 0xffffffffull);
}

TEST(log2Floor, COMMON) {
    EXPECT_EQ(0, qtr::log2Floor(1));
    EXPECT_EQ(1, qtr::log2Floor(2));
    EXPECT_EQ(1, qtr::log2Floor(3));
    EXPECT_EQ(2, qtr::log2Floor(4));
    EXPECT_EQ(2, qtr::log2Floor(5));
    EXPECT_EQ(2, qtr::log2Floor(6));
    EXPECT_EQ(2, qtr::log2Floor(7));
    EXPECT_EQ(3, qtr::log2Floor(8));
    EXPECT_EQ(10, qtr::log2Floor(1025));
}