#include <gtest/gtest.h>

#include "Fingerprint.h"

template<std::size_t FINGERPRINT_SIZE>
void testNoAlign(std::size_t expectedCountOfBlocks) {
    Fingerprint<FINGERPRINT_SIZE> fingerprint;
    EXPECT_EQ(sizeof(fingerprint), expectedCountOfBlocks * sizeof(typename Fingerprint<FINGERPRINT_SIZE>::BlockType));
}

TEST(Fingerprint, NO_ALIGN) {
    testNoAlign<0>(0);
    testNoAlign<1>(1);
    testNoAlign<2>(1);
    testNoAlign<15>(1);
    testNoAlign<16>(1);
    testNoAlign<32>(1);
    testNoAlign<33>(1);
    testNoAlign<63>(1);
    testNoAlign<64>(1);
    testNoAlign<65>(2);
    testNoAlign<123>(2);
    testNoAlign<128>(2);
    testNoAlign<129>(3);
    testNoAlign<1024>(16);
    testNoAlign<10000>(157);
}
