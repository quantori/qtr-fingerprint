#include <vector>

#include "gtest/gtest.h"

#include <cmath>
#include <stdexcept>
#include <random>
#include <set>

#include "columns_choosing/choice_functions/PearsonCorrelationChoiceFuncUtils.h"
#include "columns_choosing/choice_functions/PearsonCorrelationChoiceFunc.h"

using namespace qtr;

double pearsonCoefficient(const std::vector<bool> &x, const std::vector<bool> &y) {
    assert(x.size() == y.size());
    uint32_t xSum = 0, ySum = 0;
    for (size_t i = 0; i < x.size(); i++) {
        xSum += x[i];
        ySum += y[i];
    }
    double xMean = (double) xSum / (double) x.size();
    double yMean = (double) ySum / (double) y.size();

    if (xSum == x.size() || xSum == 0 || ySum == y.size() || ySum == 0)
        throw std::invalid_argument("x and y can not be constant columns");

    double numerator = 0;
    double denominator1 = 0, denominator2 = 0;
    for (size_t i = 0; i < x.size(); i++) {
        numerator += ((double) x[i] - xMean) * ((double) y[i] - yMean);
        denominator1 += ((double) x[i] - xMean) * ((double) x[i] - xMean);
        denominator2 += ((double) y[i] - yMean) * ((double) y[i] - yMean);
    }
    return numerator / sqrt(denominator1 * denominator2);
}

class PearsonCorrelationChoiceFuncTests : public ::testing::Test {
protected:
    std::vector<std::vector<bool>> columns;
    std::vector<bool> x, y, z, notX, trueValues, falseValues;
    std::vector<bool> a, b, notA;

    void SetUp() override {
        a = {true, false, true};
        b = {true, false, false};
        x = {true, false, true, true, false};
        y = {true, false, false, true, false};
        z = {false, true, false, false, false};
        notX = {false, true, false, false, true};
        trueValues = {true, true, true, true, true};
        falseValues = {false, false, false, false, false};
        columns = {x,
                   y,
                   trueValues,
                   notX,
                   falseValues,
                   falseValues,
                   trueValues,
                   z};
    }
};

TEST_F(PearsonCorrelationChoiceFuncTests, isConstColumnTests) {
    EXPECT_TRUE(isConstColumn(trueValues));
    EXPECT_TRUE(isConstColumn(falseValues));
    EXPECT_FALSE(isConstColumn(x));
    EXPECT_FALSE(isConstColumn(y));
    EXPECT_FALSE(isConstColumn(z));
}

TEST_F(PearsonCorrelationChoiceFuncTests, FormulaTests) {
    EXPECT_DOUBLE_EQ(pearsonCoefficient(a, b), 0.5);
    EXPECT_DOUBLE_EQ(pearsonCoefficient(x, x), 1);
    EXPECT_DOUBLE_EQ(pearsonCoefficient(x, notX), -1);
    EXPECT_DOUBLE_EQ(pearsonCoefficient(x, y), 2.0 / 3);
    EXPECT_DOUBLE_EQ(pearsonCoefficient(y, notX), -2.0 / 3);

    EXPECT_THROW(pearsonCoefficient(trueValues, x), std::invalid_argument);
    EXPECT_THROW(pearsonCoefficient(x, falseValues), std::invalid_argument);
    EXPECT_THROW(pearsonCoefficient(trueValues, trueValues), std::invalid_argument);
    EXPECT_THROW(pearsonCoefficient(falseValues, falseValues), std::invalid_argument);
}

TEST_F(PearsonCorrelationChoiceFuncTests, OptimizedFormulaTests) {
    EXPECT_DOUBLE_EQ(findPearsonCorrelation(a, b), 0.5);
    EXPECT_DOUBLE_EQ(findPearsonCorrelation(x, x), 1);
    EXPECT_DOUBLE_EQ(findPearsonCorrelation(x, notX), -1);
    EXPECT_DOUBLE_EQ(findPearsonCorrelation(x, y), 2.0 / 3);
    EXPECT_DOUBLE_EQ(findPearsonCorrelation(y, notX), -2.0 / 3);

    EXPECT_THROW(findPearsonCorrelation(trueValues, x), std::invalid_argument);
    EXPECT_THROW(findPearsonCorrelation(x, falseValues), std::invalid_argument);
    EXPECT_THROW(findPearsonCorrelation(trueValues, trueValues), std::invalid_argument);
    EXPECT_THROW(findPearsonCorrelation(falseValues, falseValues), std::invalid_argument);
}

TEST_F(PearsonCorrelationChoiceFuncTests, MaxCorrelationTest) {
    auto actual = findMaxAbsPearsonCorrelation(columns);
    std::vector<double> expected = {1,
                                    2.0 / 3,
                                    pearsonCorrelationInfinityValue,
                                    1,
                                    pearsonCorrelationInfinityValue,
                                    pearsonCorrelationInfinityValue,
                                    pearsonCorrelationInfinityValue,
                                    0.61237243569579458
    };
    EXPECT_EQ(actual.size(), expected.size());
    for (size_t i = 0; i < actual.size(); i++)
        EXPECT_DOUBLE_EQ(actual[i], expected[i]);
}

TEST_F(PearsonCorrelationChoiceFuncTests, ChoiceFuncTest) {
    IndigoFingerprint x, y, z;
    y[1] = z[1] = true;
    y[2] = z[2] = true;
    x[3] = z[3] = true;
    x[4] = y[4] = true;
    x[5] = y[5] = z[5] = true;
    IndigoFingerprintTable fingerprints;
    fingerprints.emplace_back(x);
    fingerprints.emplace_back(y);
    fingerprints.emplace_back(z);
    auto columns = PearsonCorrelationChoiceFunc()(fingerprints);
    std::set<int> actualFirstColumns = {columns[0], columns[1]};
    std::set<int> expectedFirstColumns = {3, 4};
    EXPECT_EQ(expectedFirstColumns, actualFirstColumns);
    std::set<int> actualMiddleColumns = {columns[2], columns[3]};
    std::set<int> expectedMiddleColumns = {1, 2};
    EXPECT_EQ(expectedMiddleColumns, actualMiddleColumns);

    for (size_t i = 4; i < columns.size(); i++) { // check all other columns
        EXPECT_TRUE(columns[i] == 0 || columns[i] >= 5);
    }
};

