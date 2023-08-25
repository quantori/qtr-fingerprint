#include "Profiling.h"

#include <chrono>
#include <utility>

using namespace std;
using namespace qtr;


#define PROFILING_POOL_BEGIN \
    ProfilingPool &instance = getInstance();  \
    lock_guard<mutex> lockGuard(instance._lock)

ProfilingPool &ProfilingPool::getInstance() {
    static ProfilingPool _instance;
    return _instance;
}

std::unordered_map<std::string, float> ProfilingPool::getStatistics() {
    PROFILING_POOL_BEGIN;
    return instance._statistics;
}

void ProfilingPool::addRecord(const std::string &label, float result) {
    PROFILING_POOL_BEGIN;
    instance._statistics[label] += result;
}

ProfilingTimer::ProfilingTimer(std::string label) : _label(std::move(label)), _duration(-1), _stopped(false) {
    _start = profilingClock::now();
}

ProfilingTimer::~ProfilingTimer() {
    stop();
}

float ProfilingTimer::stop() {
    if (_stopped)
        return _duration;
    _stopped = true;
    chrono::duration<float> duration = profilingClock::now() - _start;
    ProfilingPool::addRecord(_label, duration.count());
    return _duration = duration.count();
}
