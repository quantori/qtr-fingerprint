#include "Fingerprint.h"

void Fingerprint::setBytes(const std::vector<std::byte> &bytes) {
    reset();
    for(size_t i = 0; i < bytes.size(); i++)
        for(size_t j = 0; j < CHAR_BIT; j++)
            if (bool((bytes[i] >> j) & std::byte(1))) 
                set(i*CHAR_BIT + j);
}

std::vector<std::byte> Fingerprint::getBytes() const {

    std::vector<std::byte> result(fingerprintSizeInBytes, std::byte(0));
    
    for(size_t i = 0; i < result.size(); i++)
        for(size_t j = 0; j < CHAR_BIT; j++)
            if (test(i*CHAR_BIT + j)) 
                result[i] |= (std::byte(1) << j);

    return result;
}
