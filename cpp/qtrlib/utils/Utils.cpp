#include "Utils.h"

#include <stdexcept>

using namespace std;

namespace qtr {

    bool endsWith(const string &a, const string &b) {
        return a.size() >= b.size() &&
               (0 == a.compare(a.size() - b.size(), b.size(), b));
    }

    int chexToInt(char letter) {
        return letter >= '0' && letter <= '9' ? letter - '0' : letter - 'a' + 10;
    }

    vector <filesystem::path> findFiles(const filesystem::path &pathToDir, string extension) {
        vector<filesystem::path> files;
        if (!extension.empty() && extension[0] != '.')
            extension = "." + extension;
        for (const auto &entry: filesystem::recursive_directory_iterator(pathToDir)) {
            if (entry.path().stem() != ".DS_Store" && entry.path().extension() == extension) {
                files.push_back(entry.path());
            }
        }
        return files;
    }

    void askAboutContinue(const string &question) {
        cout << question << ". Continue? (y/N): ";
        string userAnswer;
        cin >> userAnswer;
        if (userAnswer != "Y" && userAnswer != "y") {
            cout << "abort\n";
            exit(-1);
        }
    }

    string generateDbName(const vector <filesystem::path> &dataDirPaths,
                          const filesystem::path &otherDataPath) {
        for (size_t i = 0;; i++) {
            bool ok = true;
            string newDbName = "db_" + to_string(i);
            for (auto &dir: dataDirPaths) {
                ok &= !filesystem::exists(dir / newDbName);
            }
            ok &= !filesystem::exists(otherDataPath / newDbName);
            if (ok)
                return newDbName;
        }
    }

    void
    initLogging(char **argv, google::LogSeverity severity, const char *base_filename, bool alsoLogToStderr) {
        google::InitGoogleLogging(argv[0]);
        google::SetLogDestination(severity, base_filename);
        FLAGS_alsologtostderr = alsoLogToStderr;
    }

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