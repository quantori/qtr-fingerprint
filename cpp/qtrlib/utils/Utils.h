#pragma once

#include <string>
#include <iostream>
#include <filesystem>
#include <vector>

#include <glog/logging.h>

namespace qtr {

/**
 * Check, if string a ends with string b.
 * @param a
 * @param b
 * @return true, if b is end of string a, false otherwise.
 */
bool endsWith(const std::string &a, const std::string &b);


/**
 * Print message if argument is empty and finish program in that case
 * @param argument
 * @param message
 */
void emptyArgument(const std::string &argument, const std::string &message);

/**
 * Return converted hex char to decimal
 */
int chexToInt(char letter);

/**
 * returns pos bit of a value, makes no checks for faster performance
 * @param value
 * @param pos
 * @return true or false
 */
inline bool getBit(uint64_t value, uint32_t pos) {
    return (value >> pos) & 1;
}

static const int BIT_IN_BYTE = 8;

/**
 * convert count of bytes to count of bits
 * @param bytes
 * @return bits
 */
constexpr inline std::size_t fromBytesToBits(std::size_t bytes) {
    return bytes * BIT_IN_BYTE;
}

/**
 * Sum of a and b should be less then std::size_t max value
 * @param a
 * @param b should be not null
 * @return Ceil of a / b
 */
constexpr inline std::size_t divideIntegersCeil(std::size_t a, std::size_t b) {
    return (a + b - 1) / b;
}

/**
 * Finds all files with extension in dir, recursively
 * @param pathToDir
 * @param extension
 * @return vector of filenames
 */
std::vector<std::filesystem::path> findFiles(const std::filesystem::path &pathToDir, const std::string &extension);

} // namespace qtr