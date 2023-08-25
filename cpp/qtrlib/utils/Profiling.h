#pragma once

#include <mutex>
#include <unordered_map>
#include "chrono"

#define CONCATENATE_TOKENS_IMPL(token1, token2) token1##token2
#define CONCATENATE_TOKENS(token1, token2) CONCATENATE_TOKENS_IMPL(token1, token2)
#define ProfileScope(label) qtr::ProfilingTimer CONCATENATE_TOKENS(profilingTimer, __LINE__) (label)

namespace qtr {

    class ProfilingPool {
    public:
        static std::unordered_map<std::string, float> getStatistics();

        static void addRecord(const std::string &label, float result);

    private:
        static ProfilingPool &getInstance();

        std::unordered_map<std::string, float> _statistics;
        std::mutex _lock;
    };

    class ProfilingTimer {
    public:
        explicit ProfilingTimer(std::string label);

        float stop();

        ~ProfilingTimer();

    private:
        using profilingClock = std::chrono::system_clock;

        std::string _label;
        decltype(profilingClock::now()) _start;
        float _duration;
        bool _stopped;
    };

} // qtr
