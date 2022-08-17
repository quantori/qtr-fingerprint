#include "gtest/gtest.h"

#include "ColumnsSelectionUtils.h"

using namespace qtr;

TEST(ColumnsChoiceUtilsTests, SortIndexByValuesTest) {
    std::map<size_t, double> values = {{0, 1.0},
                                       {1, -0.33},
                                       {2, 0.12}};
    auto actual = sortIndexesByValues(values);
    std::vector<size_t> expected = {1, 2, 0};
    EXPECT_EQ(expected, actual);
}

TEST(ColumnsChoiceUtilsTests, fingerprintsToColumnsTest) {
    IndigoFingerprint x, y, z;
    x[0] = x[1] = y[1] = z[2] = z[1] = true;
    IndigoFingerprintTable fingerprints;
    fingerprints.emplace_back(x);
    fingerprints.emplace_back(y);
    fingerprints.emplace_back(z);
    std::vector<std::vector<bool>> expected = {{true,  false, false},
                                               {true,  true,  true},
                                               {false, false, true},
                                               {false, false, false}};
    auto actualFull = fingerprintsToColumns(fingerprints);
    EXPECT_EQ(actualFull[0], expected[0]);
    EXPECT_EQ(actualFull[1], expected[1]);
    EXPECT_EQ(actualFull[2], expected[2]);
    for (size_t i = 3; i < actualFull.size(); i++) {
        EXPECT_EQ(actualFull[i], actualFull[3]);
    }

    auto actualSub = fingerprintsToColumns(fingerprints, {2, 10, 1, 0});
    EXPECT_EQ(actualSub[0], actualFull[0]);
    EXPECT_EQ(actualSub[1], actualFull[1]);
    EXPECT_EQ(actualSub[2], actualFull[2]);
    EXPECT_EQ(actualSub[10], actualFull[10]);
    EXPECT_EQ(actualSub.size(), 4);
}