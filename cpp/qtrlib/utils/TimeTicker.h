#pragma once

#include <vector>
#include <string>
#include <chrono>

namespace qtr {

    class TimeTicker {
    public:
        TimeTicker();

        double tick(const std::string &message = "");

        [[nodiscard]] double elapsedTime() const;

        void logResults() const;

    private:
        std::vector<decltype(std::chrono::high_resolution_clock::now())> _timePoints;
        std::vector <std::pair<std::string, double>> _results;
    };

} // namespace qtr
