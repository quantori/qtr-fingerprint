#include "gtest/gtest.h"
#include "../utils/TmpDirFixture.h"

#include "data_io/string_random_access_table/StringRandomAccessTable.h"
#include "string_table_io/StringTableWriter.h"

class IndexesRandomAccessTableTests : public TmpDirFixture {
public:
    std::filesystem::path getTableDir() const {
        return getTmpDir() / "table";
    }

    std::filesystem::path getTablePath() const {
        return getTmpDir() / ("stringsTable" + qtr::stringTableExtension);
    }

    void generateStringsTable() {
        qtr::StringTableWriter stWriter(getTablePath());
        for (uint64_t i = 0; i < 10; i++) {
            stWriter << std::make_pair(i, std::string(i, char('A' + i)));
        }
    }

    void SetUp() override {
        TmpDirFixture::SetUp();
        std::filesystem::create_directory(getTableDir());
        generateStringsTable();
    }
};

TEST_F(IndexesRandomAccessTableTests, AskAfterCreate) {
    qtr::StringRandomAccessTable table(getTablePath(), getTableDir());
    for (size_t i: {9, 5, 4, 1, 0, 8, 7, 6, 2, 3, 0, 9}) {
        std::string queryAnswer = table[i];
        EXPECT_EQ(i, queryAnswer.size());
        EXPECT_EQ(i, std::count(queryAnswer.begin(), queryAnswer.end(), char('A' + i)));
    }
}


TEST_F(IndexesRandomAccessTableTests, AskAfterLoad) {
    {
        qtr::StringRandomAccessTable table(getTablePath(), getTableDir());
    }
    qtr::StringRandomAccessTable table(getTableDir());
    for (size_t i: {9, 5, 4, 1, 0, 8, 7, 6, 2, 3, 0, 9}) {
        std::string queryAnswer = table[i];
        EXPECT_EQ(i, queryAnswer.size());
        EXPECT_EQ(std::string(i, char('A' + i)), queryAnswer);
    }
}

