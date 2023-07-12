#include "TimeTicker.h"

#include "glog/logging.h"

using namespace std;

namespace qtr {

    double TimeTicker::tick(const string &message) {
        _timePoints.emplace_back(chrono::high_resolution_clock::now());
        chrono::duration<double> t = _timePoints.back() - _timePoints[_timePoints.size() - 2];
        if (!message.empty()) {
            LOG(INFO) << message << ": " << t.count() << "sec";
            _results.emplace_back(message, t.count());
        }
        return t.count();
    }

    double TimeTicker::elapsedTime() const {
        return chrono::duration<double>(_timePoints.back() - _timePoints.front()).count();
    }

    void TimeTicker::logResults() const {
        for (auto &[message, duration]: _results) {
            LOG(INFO) << message << ": " << duration << "sec";
        }
        LOG(INFO) << "Elapsed time: " << elapsedTime();
    }

    TimeTicker::TimeTicker() {
        _timePoints.emplace_back(std::chrono::high_resolution_clock::now());
    }

} // namespace qtr