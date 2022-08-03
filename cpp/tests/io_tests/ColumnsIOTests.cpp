#include <filesystem>

#include "gtest/gtest.h"

#include "ColumnsIO.h"

//static const std::filesystem::path tmpColumns = std::filesystem::temp_directory_path() / "columnsTmp";
static const std::filesystem::path tmpColumns = "/home/Vsevolod.Vaskin/qtr-fingerprint/data/tests/columnsTmp";

static void writeTmpColumns(const std::vector<qtr::column_t> &values) {
    qtr::ColumnsWriter writer(tmpColumns);
    writer.write(values);
}

static std::vector<qtr::column_t> readTmpColumns() {
    qtr::ColumnsReader reader(tmpColumns);
    return reader.readAll();
}

TEST(ColumnsIOTest, EmptyFile) {
    std::vector<qtr::column_t> expected;
    writeTmpColumns(expected);
    auto actual = readTmpColumns();
    EXPECT_EQ(expected, actual);
}

TEST(ColumnsIOTest, Value) {
    std::vector<qtr::column_t> expected = {1};
    writeTmpColumns(expected);
    auto actual = readTmpColumns();
    EXPECT_EQ(expected, actual);
}

TEST(ColumnsIOTest, Values) {
    std::vector<qtr::column_t> expected = {5, 6, 4, 1, 3, 0, 2};
    writeTmpColumns(expected);
    auto actual = readTmpColumns();
    EXPECT_EQ(expected, actual);
}
