#include <fstream>
#include <numeric>

#include "gtest/gtest.h"

#include "../utils/TmpDirFixture.h"
#include "BufferedWriter.h"

class BufferedWriterTests : public TmpDirFixture {
public:
    static const size_t StrLen = 256;

    std::filesystem::path file;
    char str[StrLen];

    void SetUp() override {
        TmpDirFixture::SetUp();
        file = getTmpDir() / "file.txt";
        std::iota(str, str + StrLen, 0);
    }

    using TmpDirFixture::TearDown;

    void checkWrittenFileContainsStr() {
        std::ifstream reader(file);
        for (char i: str) {
            EXPECT_EQ((int)(unsigned char)i, (int)reader.get());
        }
        EXPECT_EQ(EOF, reader.get());
    }
};

TEST_F(BufferedWriterTests, NoWriteTest) {
    {
        qtr::BufferedWriter<10> writer(file);
    }
    std::ifstream reader(file);
    EXPECT_EQ(reader.peek(), EOF);
}
TEST_F(BufferedWriterTests, EmptyWriteTest) {
    {
        qtr::BufferedWriter<10> writer(file);
        writer.write(str, 0);
    }
    std::ifstream reader(file);
    EXPECT_EQ(reader.peek(), EOF);
}



TEST_F(BufferedWriterTests, EmptyFlushTest) {
    {
        qtr::BufferedWriter<10> writer(file);
        writer.flush();
    }
    std::ifstream reader(file);
    EXPECT_EQ(reader.peek(), EOF);
}

TEST_F(BufferedWriterTests, PutTest) {
    {
        qtr::BufferedWriter<10> writer(file);
        for (char i: str) {
            writer.put(i);
        }
    }
    checkWrittenFileContainsStr();
}

TEST_F(BufferedWriterTests, WriteTest) {
    {
        qtr::BufferedWriter<10> writer(file);
        writer.write(str, StrLen);
    }
    checkWrittenFileContainsStr();
}

TEST_F(BufferedWriterTests, BothPutAndWriteTest) {
    {
        qtr::BufferedWriter<7> writer(file);
        for (size_t i = 0; i < StrLen;) {
            if (i % 13 == 0) {
                size_t symbolsToWrite = std::min(i % 11 + 1, StrLen - i);
                writer.write(str + i, symbolsToWrite);
                i += symbolsToWrite;
            } else {
                writer.put(str[i++]);
            }
        }
    }
    checkWrittenFileContainsStr();
}





