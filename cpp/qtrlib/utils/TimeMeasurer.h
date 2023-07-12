#pragma once

#include <string>
#include <unordered_map>
#include <mutex>


namespace qtr {
    class TimeMeasurer {
    public:
        using StorageType = std::unordered_map<std::string, double>;

        StorageType::iterator begin();

        StorageType::iterator end();

        class FunctionExecutionTimer {
        public:
            FunctionExecutionTimer(TimeMeasurer &statisticCollector, std::string label);

            ~FunctionExecutionTimer();

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
}