#include <gtest/gtest.h>

#include "Fingerprint.h"
#include "FingerprintTable.h"

#include <cstring>

template<size_t FINGERPRINT_SIZE>
struct FingerprintTableTest {
    static void testRowIterator(const std::vector<std::string> &fingerprints) {
        size_t tableSize = fingerprints.size();
        FingerprintTable<FINGERPRINT_SIZE> table(tableSize);
        for (size_t i = 0; i < tableSize; ++i) {
            table[i].buildFromIndigoFingerprint(fingerprints[i].c_str());
        }

        std::size_t index = 0;
        for (const Fingerprint<FINGERPRINT_SIZE> &fingerprint: table) {
            EXPECT_TRUE(fingerprint._data == table[index]._data);
            index++;
        }
    }
};

TEST(FingerprintTableRowIteratorTest, COMMON) {
    FingerprintTableTest<16 * 4>::testRowIterator({"f79bd5571c488178", "06ebefc68f856039", "5bc8f9b5f001c7dc",
                                                   "8fd86b2b8721fb5f", "a8938805d3ab99a5", "ae75187829ac7de7"});
    FingerprintTableTest<32 * 4>::testRowIterator(
            {"e3a551013fc2c16e3e4c65329b2f163d", "d7b4b75bd79faee3ba2a0d93f6c545f8"});
    FingerprintTableTest<48 * 4>::testRowIterator({"68a1252b85d6cf6a40fea7d17856fa100cf981e7a3172af7",
                                                   "1894d61203c7273ef56465509cefbf24d325ed58a63d78cc"});
}


