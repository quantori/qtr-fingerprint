#include "Utils.h"

namespace qtr {

    bool endsWith(const std::string &a, const std::string &b) {
        return a.size() >= b.size() &&
               (0 == a.compare(a.size() - b.size(), b.size(), b));
    }

    int chexToInt(char letter) {
        return letter >= '0' && letter <= '9' ? letter - '0' : letter - 'a' + 10;
    }

    std::vector<std::filesystem::path> findFiles(const std::filesystem::path &pathToDir, const std::string &extension) {
        std::vector<std::filesystem::path> sdfFiles;
        for (const auto &entry: std::filesystem::recursive_directory_iterator(pathToDir)) {
            if (entry.path().extension() == extension) {
                //            LOG(INFO) << entry.path().string() << std::endl;
                sdfFiles.push_back(entry.path());
            }
        }
        return sdfFiles;
    }

    template<>
    void emptyArgument<uint64_t>(const uint64_t &argument, const std::string &message) {
        if (argument == 0) {
            LOG(ERROR) << message;
            exit(-1);
        }
    }

    void askAboutContinue(const std::string &question) {
        std::cout << question << ". Continue? (Y/n): ";
        std::string userAnswer;
        std::cin >> userAnswer;
        if (userAnswer != "Y" && userAnswer != "y") {
            std::cout << "abort\n";
            exit(-1);
        }
    }

    std::string generateDbName(const std::vector<std::filesystem::path> &dataDirPaths,
                               const std::filesystem::path &otherDataPath) {
        for (size_t i = 0;; i++) {
            bool ok = true;
            std::string newDbName = "db_" + std::to_string(i);
            for (auto &dir: dataDirPaths) {
                ok &= !std::filesystem::exists(dir / newDbName);
            }
            ok &= !std::filesystem::exists(otherDataPath / newDbName);
            if (ok)
                return newDbName;
        }
    }

    double TimeTicker::tick(const std::string &message) {
        _timePoints.emplace_back(std::chrono::high_resolution_clock::now());
        std::chrono::duration<double> t = _timePoints.back() - _timePoints[_timePoints.size() - 2];
        LOG(INFO) << message << ": " << t.count() << "sec";
        _results.emplace_back(message, t.count());
        return t.count();
    }

    double TimeTicker::elapsedTime() const {
        return std::chrono::duration<double>(_timePoints.back() - _timePoints.front()).count();
    }

    void TimeTicker::logResults() const {
        for (auto &[message, duration]: _results) {
            LOG(INFO) << message << ": " << duration << "sec";
        }
        LOG(INFO) << "Elapsed time: " << elapsedTime();
    }


} // namespace qtr