#include "gtest/gtest.h"

#include "HuffmanCoder.h"

class HuffmanTests : public ::testing::Test {
protected:

    void testStrings(const std::vector<std::string> &strings) {
        qtr::HuffmanCoder::Builder builder;
        for (auto &s: strings) {
            builder += s;
        }
        auto coder = builder.build();
        for (auto &s: strings) {
            auto code = coder.encode(s);
            EXPECT_LE(code.size(), s.size() * CHAR_BIT);
            auto t = coder.decode(code);
            EXPECT_EQ(s, t);
        }
    }
};

TEST_F(HuffmanTests, allSymbols) {
    std::string s;
    for (int i = 1; i < 256; i++) {
        int priority = (i * 5 + 13) % 19 + 1;
        s += std::string(priority, char(i));
    }
    testStrings({s});
}

TEST_F(HuffmanTests, randomString) {
    std::string s = "This is sample text to test huffman coder on string. 1234567890."
                    "AND WE ALSO SHOULD TEST SOME CAPS LOCK SYMBOLS";
    testStrings({s});
}

TEST_F(HuffmanTests, smilesStrings) {
    testStrings({"CNC1=SC=CN1",
                 "O=S(=O)(NC1=CC=CC=C1)C1=C2C=CN=CC2=CC=C1",
                 "C[C@](O)([C@H]1CCCC2=Cc3c(C[C@]12C)cnn3-c4ccc(F)cc4)c5ccsc5",
                 "[H]N1C([H])=C2C(=C1O)C3=C(NC4=C3C([H])=C([H])C([H])=C4[H])C5=C2C6=C(N5)C([H])=C([H])C([H])=C6[H]",
                 "CCN1C(=O)C(C)=CC2=CN=C(C)N=C12",
                 "C1CCC2C(C1)CCC3=CC=CC=C23",
                 "CC(C)(N)CC1=CC=CC=C1",
                 "CC(O)CC(C)(C)C(CCC#C)C(C)(C)C"});
}