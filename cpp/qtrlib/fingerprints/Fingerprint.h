#pragma once

#include <bitset>
#include <climits>
#include <cstddef>
#include <vector>

class QtrIndigoFingerprint;

namespace qtr {

template<size_t fingerprintSizeInBytes>
class Fingerprint : public std::bitset<CHAR_BIT*fingerprintSizeInBytes>
{
public:
    static constexpr size_t sizeInBytes = fingerprintSizeInBytes;

    Fingerprint() = default;

    explicit Fingerprint(const QtrIndigoFingerprint &f);

    /**
     * Builds fingerprint from HEX string(for example - "a4210f")
     * @param s - string with fingerprint, strlen(s) should be equals to fingerprintSizeInBytes * 2
     */
    explicit Fingerprint(const std::string &s) {
        for (size_t i = 0; i < s.size(); ++i) {
            size_t j = i * 4ull;
            int currentSym = (s[i] >= 'a' && s[i] <= 'z') ? (s[i] - 'a' + 10) : (s[i] - '0');
            this->operator[](j) = currentSym & 1;
            this->operator[](j + 1) = currentSym & 2;
            this->operator[](j + 2) = currentSym & 4;
            this->operator[](j + 3) = currentSym & 8;
        }
    }

    void setBytes(const std::vector<std::byte> &bytes) {
        this->reset();
        for(size_t i = 0; i < bytes.size(); i++)
            for(size_t j = 0; j < CHAR_BIT; j++)
                if (bool((bytes[i] >> j) & std::byte(1))) 
                    this->set(i*CHAR_BIT + j);
    }

    std::vector<std::byte> getBytes() const {
        std::vector<std::byte> result(fingerprintSizeInBytes, std::byte(0));

        for(size_t i = 0; i < result.size(); i++)
            for(size_t j = 0; j < CHAR_BIT; j++)
                if (this->test(i*CHAR_BIT + j)) 
                    result[i] |= (std::byte(1) << j);

        return result;
    }

    void saveBytes(std::ostream& out) const {
        auto toWrite = getBytes();
        out.write((char*)toWrite.data(), toWrite.size());
    }
};

using IndigoFingerprint = Fingerprint<323>;
using FullIndigoFingerprint = Fingerprint<467>;

} // namespace qtr