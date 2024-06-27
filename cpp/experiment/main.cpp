#include <glog/logging.h>

#include <fstream>
#include <future>
#include <thread>

#include "ExperimentArgs.h"
#include "SearchEngineFactory.h"
#include "ExperimentStat.h"

#include "GraphMol/GraphMol.h"
#include "GraphMol/SmilesParse/SmilesParse.h"


std::vector<std::string> parseMoleculesFromFile(const std::filesystem::path &moleculesFile) {
    std::ifstream in(moleculesFile);
    std::vector<std::string> result;
    while (in.peek() != EOF) {
        std::string smiles;
        in >> smiles;
        std::string lineEnding;
        std::getline(in, lineEnding);
        if (smiles.empty()) {
            LOG(ERROR) << "Wrong formatting in queries file";
        }
        result.push_back(smiles);
    }
    LOG(INFO) << "Loaded " << result.size() << " molecules from file " << moleculesFile;
    return result;
}

struct ExperimentInfo {
    bool stopFlag;
    std::atomic_flag experimentFinished;
    double timeLimit;
    std::chrono::time_point<std::chrono::system_clock, std::chrono::system_clock::duration> startTime;
    std::mutex mutex;

    void checkTimeout() {
        std::lock_guard<std::mutex> lockGuard(mutex);
        if (stopFlag) {
            return;
        }
        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = now - startTime;
        bool timeout = duration.count() > timeLimit;
        stopFlag = timeout;
        if (timeout) {
            LOG(INFO) << "Timout detected";
        }
    }

    void startExperiment() {
        std::lock_guard<std::mutex> lockGuard(mutex);
        stopFlag = false;
        startTime = std::chrono::high_resolution_clock::now();
    }

    [[nodiscard]] double finishExperiment() const {
        std::chrono::duration<double> duration = std::chrono::high_resolution_clock::now() - startTime;
        return duration.count();
    }
};

void checkTimeout(ExperimentInfo &info) {
    // Duration is bounded by 0.1 and 1.0 to prevent too often and too rare checks
    std::chrono::duration<double> sleepDuration(std::max(0.1, std::min(1.0, info.timeLimit)));
    while (!info.experimentFinished.test()) {
        std::this_thread::sleep_for(sleepDuration);
        info.checkTimeout();
    }
}

void conductExperiment(SearchEngine &searchEngine,
                       const std::vector<std::string> &queries,
                       int maxResults,
                       double timeLimit,
                       std::ofstream &statOut) {
    ExperimentInfo info = {
            .stopFlag = false,
            .experimentFinished = false,
            .timeLimit = timeLimit,
    };
    auto checkTimeoutThread = std::thread(checkTimeout, std::ref(info));
    ExperimentStat stat;
    for (auto &query: queries) {
        LOG(INFO) << "Start " << query << " processing";
        info.startExperiment();
        auto matches = searchEngine.getMatches(query, maxResults, info.stopFlag);

        auto experimentDuration = info.finishExperiment();
        stat.add(ExperimentStat::Entity{
                .duration = experimentDuration,
                .resultsFound = (int) matches.size()
        });
        LOG(INFO) << "Finished " << query << " processing."
                  << "\n\tDuration: " << experimentDuration
                  << "\n\tresultsFound: " << matches.size();
    }
    info.experimentFinished.test_and_set();
    LOG(INFO) << "Wait timeout checker to finish";
    checkTimeoutThread.join();
    statOut << stat << std::endl;
}

int main(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    google::SetLogDestination(google::INFO, "experiment.info");
    FLAGS_alsologtostderr = true;

    ExperimentArgs args(argc, argv);

    auto searchEngine = SearchEngineFactory::create(args.searchEngineType, args.datasetDir);
    auto queries = parseMoleculesFromFile(args.queriesFile);
    std::ofstream statOut(args.statisticsFile);
    conductExperiment(*searchEngine, queries, args.maxResults, args.timeLimit, statOut);
    return 0;
}