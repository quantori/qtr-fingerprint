#pragma once

#include <string>
#include <iostream>
#include <filesystem>
#include <vector>
#include <unordered_map>

#include <glog/logging.h>
#include "glog/log_severity.h"

namespace qtr {

    class TimeTicker {
    public:
        TimeTicker() {
            _timePoints.emplace_back(std::chrono::high_resolution_clock::now());
        }

        double tick(const std::string &message = "");

        [[nodiscard]] double elapsedTime() const;

        void logResults() const;

    private:
        std::vector<decltype(std::chrono::high_resolution_clock::now())> _timePoints;
        std::vector<std::pair<std::string, double>> _results;
    };

    class TimeMeasurer {
    public:
        using StorageType = std::unordered_map<std::string, double>;

        StorageType::iterator begin();

        StorageType::iterator end();

        class FunctionTimeMeasurer {
        public:
            FunctionTimeMeasurer(TimeMeasurer &statisticCollector, std::string label);

            ~FunctionTimeMeasurer();

        private:
            TimeMeasurer &_statisticCollector;
            std::string _label;
        };

        void start(const std::string &label);

        void finish(const std::string &label);

    private:
        std::unordered_map<std::string, double> _measurements;
        std::unordered_map<std::string, decltype(std::chrono::high_resolution_clock::now())> _startPoints;
        std::mutex _lock;
    };

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
    template<typename T>
    void checkEmptyArgument(const T &argument, const std::string &message) {
        if (argument.empty()) {
            LOG(ERROR) << message;
            exit(-1);
        }
    }

    std::string generateDbName(const std::vector<std::filesystem::path> &dataDirPaths,
                               const std::filesystem::path &otherDataPath);

    void askAboutContinue(const std::string &question);

    template<>
    void checkEmptyArgument<uint64_t>(const uint64_t &argument, const std::string &message);

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
    findFiles(const std::filesystem::path &pathToDir, std::string extension = "");

    /**
     * Initialize google logging
     */
    void initLogging(char **argv, google::LogSeverity severity, const char *base_filename, bool alsoLogToStderr);
} // namespace qtr