#include "gtest/gtest.h"
#include "../utils/TmpDirFixture.h"

#include "smiles_table_io/SmilesRandomAccessTable.h"
#include "smiles_table_io/SmilesTableWriter.h"

class IndexesRandomAccessTableTests : public TmpDirFixture {
public:
    std::filesystem::path getTableDir() const {
        return getTmpDir() / "table";
    }

    std::filesystem::path getSmilesTablePath() const {
        return getTmpDir() / ("_smilesTable" + qtr::smilesTableExtension);
    }

    void generateSmilesTable() {
        qtr::SmilesTableWriter stWriter(getSmilesTablePath());
        for (uint64_t i = 0; i < 10; i++) {
            stWriter << std::make_pair(i, std::string(i, char('A' + i)));
        }
    }

    void SetUp() override {
        TmpDirFixture::SetUp();
        std::filesystem::create_directory(getTableDir());
        generateSmilesTable();
    }
};

TEST_F(IndexesRandomAccessTableTests, AskAfterCreate) {
    qtr::SmilesRandomAccessTable table(getSmilesTablePath(), getTableDir());
    for (size_t i: {9, 5, 4, 1, 0, 8, 7, 6, 2, 3, 0, 9}) {
        std::string queryAnswer = table[i];
        EXPECT_EQ(i, queryAnswer.size());
        EXPECT_EQ(i, std::count(queryAnswer.begin(), queryAnswer.end(), char('A' + i)));
    }
}


TEST_F(IndexesRandomAccessTableTests, AskAfterLoad) {
    {
        qtr::SmilesRandomAccessTable table(getSmilesTablePath(), getTableDir());
    }
    qtr::SmilesRandomAccessTable table(getTableDir());
    for (size_t i: {9, 5, 4, 1, 0, 8, 7, 6, 2, 3, 0, 9}) {
        std::string queryAnswer = table[i];
        EXPECT_EQ(i, queryAnswer.size());
        EXPECT_EQ(std::string(i, char('A' + i)), queryAnswer);
    }
}

