#include "gtest/gtest.h"

#include "../utils/TmpDirFixture.h"
#include "MemoryMapFile.h"

class MemoryMapFileTests : public TmpDirFixture {

};

TEST_F(MemoryMapFileTests, MemoryMapFileTest) {
    std::filesystem::path testFile = getTmpDir() / "mmf";
    const int expectedCapacity = 1025;
    {
        auto mmf = qtr::MemoryMapFile::create(testFile, expectedCapacity);
        EXPECT_EQ(mmf.capacity(), expectedCapacity);
        EXPECT_EQ(mmf.size(), 0);
        for (int i = 0; i < 256; i++) {
            uint32_t *value = static_cast<uint32_t *>(mmf.allocate(4));
            ASSERT_NE(value, nullptr);
            EXPECT_EQ(mmf.capacity(), expectedCapacity);
            EXPECT_EQ(mmf.size(), (i + 1) * 4);
            *value = i;
        }
        EXPECT_EQ(mmf.allocate(2), nullptr);
    }
    {
        auto mmf = qtr::MemoryMapFile::open(testFile);
        EXPECT_EQ(mmf.capacity(), expectedCapacity);
        EXPECT_EQ(mmf.size(), expectedCapacity - 1);
        char *symbol = static_cast<char *>(mmf.allocate(1));
        ASSERT_NE(symbol, nullptr);
        *symbol = 'b';
        EXPECT_EQ(mmf.allocate(1), nullptr);
    }
    {
        auto mmf = qtr::MemoryMapFile::open(testFile);
        char *symbol = static_cast<char *>(mmf.ptr(1024));
        ASSERT_NE(symbol, nullptr);
        EXPECT_EQ(*symbol, 'b');
    }
}