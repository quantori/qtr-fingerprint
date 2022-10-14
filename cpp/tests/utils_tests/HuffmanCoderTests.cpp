#include "gtest/gtest.h"
#include "../utils/TmpDirFixture.h"

#include "HuffmanCoder.h"
#include <fstream>

class HuffmanTests : public TmpDirFixture {
protected:
    std::filesystem::path codeFile;

    void SetUp() override {
        TmpDirFixture::SetUp();
        codeFile = getTmpDir() / "initFile";
    }

    void analyzeString(const std::string &s) {
        std::map<char, uint64_t> frequency;
        for (char symbol: s + " ") {
            frequency[symbol]++;
        }
        std::ofstream out(codeFile);
        for (auto &[symbol, rate]: frequency) {
            out << (int) (uint8_t) symbol << ' ' << rate << "\n";
        }
    }

    void testString(const std::string &s) {
        analyzeString(s);
        qtr::HuffmanCoder coder(codeFile);
        auto code = coder.encode(s);
        EXPECT_LE(code.size(), s.size() * CHAR_BIT);
        auto t = coder.decode(code);
        EXPECT_EQ(s, t);
    }
};

TEST_F(HuffmanTests, emptyString) {
    testString("");
}

TEST_F(HuffmanTests, allSymbols) {
    std::string s;
    for (int i = 1; i < 256; i++) {
        s += char(i);
    }
    testString(s);
}

TEST_F(HuffmanTests, randomString) {
    std::string s = "This is some text to test huffman coder on random string. 1234567890."
                    "AND WE ALSO SHOULD TEST SOME CAPS LOCK SYMBOLS";
    testString(s);
}