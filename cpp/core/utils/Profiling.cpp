#include "Profiling.h"

#include <utility>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;

#define PROFILING_POOL_BEGIN \
    ProfilingPool &instance = getInstance();  \
    lock_guard<mutex> lockGuard(instance._lock)

ProfilingPool &ProfilingPool::getInstance() {
    static ProfilingPool _instance;
    return _instance;
}

unordered_map<string, float> ProfilingPool::getStatistics() {
    PROFILING_POOL_BEGIN;
    return instance._statistics;
}

void ProfilingPool::addRecord(const string &label, float result) {
    auto start = std::chrono::system_clock::now();
    PROFILING_POOL_BEGIN;
    instance._statistics[label] += result;
    chrono::duration<float> duration = chrono::system_clock::now() - start;
    instance._addRecordTimer += duration.count();
}

void ProfilingPool::showStatistics(ostream &out) {
    PROFILING_POOL_BEGIN;
    auto statistics = instance._statistics;
    out << "Profiling information:\n";
    vector<pair<string, float>> stat(statistics.begin(), statistics.end());
    stat.emplace_back("Time spent on profiling: ", instance._addRecordTimer);
    sort(stat.begin(), stat.end());
    for (auto &[label, time]: stat) {
        out << fixed << setprecision(6) << label << ": " << time << " seconds\n";
    }
}

ProfilingTimer::ProfilingTimer(string label) : _label(std::move(label)), _duration(-1), _stopped(false) {
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