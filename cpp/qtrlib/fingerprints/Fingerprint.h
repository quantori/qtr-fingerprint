#pragma once

#include <bitset>
#include <climits>
#include <cstddef>
#include <vector>

constexpr size_t fingerprintSizeInBytes = 467;

class Fingerprint : public std::bitset<CHAR_BIT*fingerprintSizeInBytes>
{
public:
    void setBytes(const std::vector<std::byte> &bytes);
    std::vector<std::byte> getBytes() const;
};
