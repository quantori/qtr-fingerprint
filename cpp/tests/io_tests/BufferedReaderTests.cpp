#include <fstream>

#include "gtest/gtest.h"

#include "../utils/TmpDirFixture.h"
#include "BufferedReader.h"

class BufferedReaderTests : public TmpDirFixture {
public:
    static const size_t StrLen = 256;

    std::filesystem::path fileQtr, fileStd;
    std::string str;

    void SetUp() override {
        TmpDirFixture::SetUp();
        fileQtr = getTmpDir() / "fileQtr.txt";
        fileStd = getTmpDir() / "fileStd.txt";
        for (int i = 0; i < StrLen; i++)
            str += char(i);
    }

    void writeEmpty() {
        std::ofstream writerStd(fileStd), writerQtr(fileQtr);
    }

    void writeStr() {
        std::ofstream writerStd(fileStd), writerQtr(fileQtr);
        writerStd << str;
        writerQtr << str;
    }
};

TEST_F(BufferedReaderTests, EmptyFileEOFTests) {
    writeEmpty();
    std::ifstream readerStd(fileStd);
    qtr::BufferedReader<10> readerQrt(fileQtr);
    EXPECT_EQ(readerQrt.eof(), readerStd.eof());
    readerQrt.get();
    readerStd.get();
    EXPECT_EQ(readerQrt.eof(), readerStd.eof());
}

TEST_F(BufferedReaderTests, GetTest) {
    writeStr();
    qtr::BufferedReader<10> reader(fileQtr);
    for (char c: str) {
        EXPECT_EQ((int) (unsigned char) c, reader.get());
    }
    EXPECT_FALSE(reader.eof());
    EXPECT_EQ(EOF, reader.get());
    EXPECT_TRUE(reader.eof());
}

TEST_F(BufferedReaderTests, ReadAllTest) {
    writeStr();
    qtr::BufferedReader<10> reader(fileQtr);
    char buf[StrLen];
    reader.read(buf, StrLen);
    EXPECT_EQ(0, std::memcmp(buf, str.c_str(), StrLen));
    EXPECT_FALSE(reader.eof());
    EXPECT_EQ(EOF, reader.get());
    EXPECT_TRUE(reader.eof());
}

TEST_F(BufferedReaderTests, ReadMoreThanAllTest) {
    writeStr();
    qtr::BufferedReader<10> reader(fileQtr);
    char buf[2 * StrLen];
    std::fill(buf, buf + 2 * StrLen, char(255));
    reader.read(buf, StrLen);
    EXPECT_EQ(0, std::memcmp(buf, str.c_str(), StrLen));
    size_t strLen = StrLen; // StrLen doesn't work in #define EXPECT_EQ (undefined reference), so we create a copy of it
    EXPECT_EQ(strLen, std::count(buf + strLen, buf + strLen * 2, char(255)));
    EXPECT_FALSE(reader.eof());
    EXPECT_EQ(EOF, reader.get());
    EXPECT_TRUE(reader.eof());
}

TEST_F(BufferedReaderTests, ReadByBlocksTest) {
    writeStr();
    qtr::BufferedReader<10> reader(fileQtr);
    char buf[2 * StrLen];
    std::fill(buf, buf + 2 * StrLen, char(255));
    for (size_t i = 0; i < StrLen; i += 7) {
        reader.read(buf + i, 7);
    }
    size_t strLen = StrLen; // StrLen doesn't work in #define EXPECT_EQ (undefined reference), so we create a copy of it
    EXPECT_EQ(0, std::memcmp(buf, str.c_str(), strLen));
    EXPECT_EQ(strLen, std::count(buf + strLen, buf + strLen * 2, char(255)));
    EXPECT_TRUE(reader.eof());
}

TEST_F(BufferedReaderTests, ReadAndGetTest) {
    writeStr();
    qtr::BufferedReader<10> reader(fileQtr);
    char buf[7];
    for (size_t i = 0; i < StrLen;) {
        if (i % 11 == 0) {
            size_t readCnt = i % 7 + 1;
            reader.read(buf, readCnt);
            EXPECT_EQ(0, std::memcmp(buf, str.c_str() + i, readCnt));
            i += readCnt;
        } else {
            EXPECT_EQ((int) (unsigned char) str[i++], reader.get());
        }
    }
    EXPECT_FALSE(reader.eof());
    EXPECT_EQ(EOF, reader.get());
    EXPECT_TRUE(reader.eof());
}