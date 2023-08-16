#include "TimeMeasurer.h"

#include "glog/logging.h"

using namespace std;

namespace qtr {
    TimeMeasurer::StorageType::iterator TimeMeasurer::begin() {
        return _measurements.begin();
    }

    unordered_map<string, double>::iterator TimeMeasurer::end() {
        return _measurements.end();
    }

    void TimeMeasurer::start(const string &label) {
        LOG(INFO) << "Start " << label;
        lock_guard<mutex> lock(_lock);
        if (_measurements.contains(label) || _startPoints.contains(label))
            throw std::invalid_argument("Such label already exists: " + label);
        _startPoints[label] = chrono::high_resolution_clock::now();
    }

    void TimeMeasurer::finish(const string &label) {
        lock_guard<mutex> lock(_lock);
        if (!_startPoints.contains(label))
            throw std::invalid_argument("Such label does not exist: " + label);
        chrono::duration<double> duration = chrono::high_resolution_clock::now() - _startPoints[label];
        _measurements[label] = duration.count();
        LOG(INFO) << "Finish " << label;
    }

    TimeMeasurer::FunctionExecutionTimer::FunctionExecutionTimer(TimeMeasurer &statisticCollector,
                                                                 std::string label)
            : _statisticCollector(statisticCollector), _label(std::move(label)) {
        _statisticCollector.start(_label);
    }

    TimeMeasurer::FunctionExecutionTimer::~FunctionExecutionTimer() {
        _statisticCollector.finish(_label);
    }
} // namespace qtr
