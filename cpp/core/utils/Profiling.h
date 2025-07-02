#pragma once

#include <mutex>
#include <unordered_map>
#include <filesystem>
#include <chrono>
#include <ostream>
#include <string>


#ifdef QTR_PROFILING
#define CONCATENATE_TOKENS_IMPL(token1, token2) token1##token2
#define CONCATENATE_TOKENS(token1, token2) CONCATENATE_TOKENS_IMPL(token1, token2)

#define CreateProfiler(timerName, label) ProfilingTimer timerName("[" + std::filesystem::path(__FILE__).filename().string() + ":" + std::to_string(__LINE__) + "] " + (label))
#define StopProfiler(timerName) (timerName).stop()
#define ProfileScope(label) CreateProfiler(CONCATENATE_TOKENS(profilingTimer, __LINE__), label)

#else
#define CreateProfiler(profilerName, label)
#define ProfileScope(label)
#define StopProfiler(profilerName)
#endif

class ProfilingPool {
public:
    static std::unordered_map<std::string, float> getStatistics();

    static void showStatistics(std::ostream &out);

    static void addRecord(const std::string &label, float result);

private:
    static ProfilingPool &getInstance();

    std::unordered_map<std::string, float> _statistics;
    std::mutex _lock;
    float _addRecordTimer = 0;
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
