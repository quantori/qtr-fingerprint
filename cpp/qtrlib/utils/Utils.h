#pragma once

#include <string>
#include <iostream>
#include <filesystem>
#include <vector>
#include <unordered_map>
#include <type_traits>

#include <glog/logging.h>
#include "glog/log_severity.h"

namespace qtr {
    /**
     * Check, if string a ends with string b.
     * @param a
     * @param b
     * @return true, if b is end of string a, false otherwise.
     */
    bool endsWith(const std::string &a, const std::string &b);


    template<typename T>
    class HasEmptyMethod {
        typedef char hasEmpty;
        struct hasNotEmpty {
            char x[2];
        };

        template<typename C>
        static hasEmpty test(decltype(&C::empty)) { return hasEmpty(); }

        template<typename C>
        static hasNotEmpty test(...) { return hasNotEmpty(); }

    public:
        enum {
            value = sizeof(test<T>(0)) == sizeof(char)
        };
    };

    template<typename T, bool hasEmpty, typename D>
    struct EmptyChecker {
        D _emptyVal;

        explicit EmptyChecker(const D &emptyVal) : _emptyVal(emptyVal) {};

        bool check(const T &val) {
            return val == _emptyVal;
        }
    };

    template<typename T, typename D>
    struct EmptyChecker<T, true, D> {

        explicit EmptyChecker(const D &) {};

        static bool check(const T &val) {
            return val.empty();
        }
    };

    template<typename T, typename D = T>
    bool checkEmpty(const T &val, const D &emptyVal) {
        return EmptyChecker<T, HasEmptyMethod<T>::value, D>(emptyVal).check(val);
    }

    /**
     * Print message if argument is empty and finish program in that case
     * @param argument
     * @param message
     */
    template<typename T>
    void checkEmptyArgument(const T &argument, const std::string &message) {
        if (checkEmpty<T>(argument, T())) {
            LOG(ERROR) << message;
            exit(-1);
        }
    }

    void askAboutContinue(const std::string &question);

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

    template<typename T>
    constexpr inline T lowerOrderBits(T number, size_t bits_count) {
        return number & ((T(1) << bits_count) - 1);
    }

    constexpr inline size_t log2Floor(size_t number) {
        size_t answer = 0;
        while (number > 1) {
            answer++;
            number >>= 1u;
        }
        return answer;
    }

    /**
     * Finds all files with extension in dir, recursively
     * @param pathToDir
     * @param extension
     * @return vector of filenames
     */
    std::vector<std::filesystem::path>
    findFiles(const std::filesystem::path &pathToDir, std::string extension);

    /**
     * Initialize google logging
     */
    void initLogging(char **argv, google::LogSeverity severity, const char *base_filename, bool alsoLogToStderr);

    /**
     * @param lhs
     * @param rhs
     * @return concatenation of @a lhs and @a rhs
     */
    template<typename T, typename D>
    std::vector<T> concatVectors(const std::vector<T> &lhs, const std::vector<D> &rhs) {
        std::vector<T> res = lhs;
        std::copy(rhs.begin(), rhs.end(), std::back_inserter(res));
        return res;
    }

    template<typename Enum>
    auto makeStringToEnumFunction(const std::unordered_map<std::string, Enum> &mapping, Enum badValue) {
        return [&mapping, badValue](const std::string &s) {
            auto it = mapping.find(s);
            if (it == mapping.end())
                return badValue;
            return it->second;
        };
    }

    template<typename T, typename = void>
    struct Printable : std::false_type {
    };

    template<typename T>
    struct Printable<T, std::void_t<decltype(std::cout << std::declval<T>())>> : std::true_type {
    };

    template<typename T>
    std::string toString(const T &value) {
        std::ostringstream oss;
        if constexpr (Printable<T>::value) {
            oss << value;
        } else {
            oss << "[";
            for (const auto &element: value) {
                oss << toString(element) << ", ";
            }
            std::string result = oss.str();
            if (result.size() > 1) {
                result.pop_back();
                result.pop_back();
            }
            result.push_back(']');
            return result;
        }
        return oss.str();
    }

    [[noreturn]] void logErrorAndExit(const std::string &message);
} // namespace qtr