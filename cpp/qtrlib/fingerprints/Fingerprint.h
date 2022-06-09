#pragma once

#include <bitset>
#include <climits>
#include <cstddef>
#include <vector>

template<size_t fingerprintSizeInBytes>
class Fingerprint : public std::bitset<CHAR_BIT*fingerprintSizeInBytes>
{
public:
    static constexpr size_t sizeInBytes = fingerprintSizeInBytes;

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
};

using FingerprintForIndigo = Fingerprint<467>;
