#include <gtest/gtest.h>

#include "Utils.h"

#include <bitset>

using namespace std;

static const int VALUE_SIZE_BITS = fromBytesToBits(sizeof(uint64_t));

TEST(getBit, ALL_ZEROES) {
    bitset<VALUE_SIZE_BITS> curr = 0;
    for (int i = 0; i < VALUE_SIZE_BITS; ++i) {
        EXPECT_EQ(curr[i], 0);
        EXPECT_EQ(curr[i], getBit(curr.to_ullong(), i));
    }
}

TEST(getBit, ALL_ONES) {
    bitset<VALUE_SIZE_BITS> curr = -1;
    for (int i = 0; i < VALUE_SIZE_BITS; ++i) {
        EXPECT_EQ(curr[i], 1);
        EXPECT_EQ(curr[i], getBit(curr.to_ullong(), i));
    }
}

TEST(getBit, RANDOM) {
    bitset<VALUE_SIZE_BITS> curr = 17432891238;
    for (int i = 0; i < VALUE_SIZE_BITS; ++i) {
        EXPECT_EQ(curr[i], getBit(curr.to_ullong(), i));
    }
}

TEST(fromBytesToBits, COMMON) {
    EXPECT_EQ(fromBytesToBits(0), 0);
    EXPECT_EQ(fromBytesToBits(1), 8);
    EXPECT_EQ(fromBytesToBits(1000), 8000);
}

TEST(divideIntegersCeil, COMMON) {
    EXPECT_EQ(divideIntegersCeil(0, 5), 0);
    EXPECT_EQ(divideIntegersCeil(1, 2), 1);
    EXPECT_EQ(divideIntegersCeil(1, 3), 1);
    EXPECT_EQ(divideIntegersCeil(3, 3), 1);
    EXPECT_EQ(divideIntegersCeil(4, 3), 2);
    EXPECT_EQ(divideIntegersCeil(15, 4), 4);
}