#include "gtest/gtest.h"

#include "ColumnsChoiceUtils.h"

using namespace qtr;

TEST(ColumnsChoiceUtilsTests, SortIndexByValuesTest) {
    std::vector<double> values = {1.0, -0.33, 0.12};
    auto actual = sortIndexesByValues(values);
    std::vector<int> expected = {1, 2, 0};
    EXPECT_EQ(expected, actual);
}

TEST(ColumnsChoiceUtilsTests, fingerprintsToColumnsTest) {
    IndigoFingerprint x, y, z;
    x[0] = x[1] = y[1] = z[2] = z[1] = true;
    IndigoFingerprintTable fingerprints;
    fingerprints.emplace_back(x);
    fingerprints.emplace_back(y);
    fingerprints.emplace_back(z);
    auto columns = fingerprintsToColumns(fingerprints);
    std::vector<std::vector<bool>> expected = {{true,  false, false},
                                               {true,  true,  true},
                                               {false, false, true},
                                               {false, false, false}};
    EXPECT_EQ(columns[0], expected[0]);
    EXPECT_EQ(columns[1], expected[1]);
    EXPECT_EQ(columns[2], expected[2]);
    for (size_t i = 3; i < columns.size(); i++) {
        EXPECT_EQ(columns[i], columns[3]);
    }
}