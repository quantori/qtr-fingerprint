#include <glog/logging.h>

#include <fstream>
#include <future>
#include <thread>
#include <execution>

#include "ExperimentArgs.h"
#include "ExperimentStat.h"
#include "SearchEngineConcept.h"
#include "RDKitSearchEngine.h"
#include "IndigoSearchEngine.h"
#include "BallTreeSearchEngine.h"
#include "TrieSearchEngine.h"
#include "IndigoBruteForceSearchEngine.h"
#include "QueriesParser.h"
#include "SmilesDirParser.h"

#include "GraphMol/GraphMol.h"
#include "GraphMol/SmilesParse/SmilesParse.h"


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
    // Duration is bounded by 0.001 and 0.01 to prevent too often and too rare checks
    std::chrono::duration<double> sleepDuration(std::max(0.001, std::min(0.01, info.timeLimit)));
    while (!info.experimentFinished.test()) {
        std::this_thread::sleep_for(sleepDuration);
        info.checkTimeout();
    }
}

template<SearchEngine SE>
void conductExperiment(SE &searchEngine,
                       const std::vector<std::string> &queries,
                       int maxResults,
                       double timeLimit,
                       std::ofstream &statOut) {
    ExperimentStat stat;
    std::mutex statMutex;
    std::atomic_int64_t counter = 0;
    std::for_each(std::execution::par, queries.begin(), queries.end(), [&](const std::string &query) {
        auto i = counter++;
        LOG(INFO) << "Start " << query << " processing (" << i + 1 << ")";
        decltype(searchEngine.smilesToQueryMolecule(query)) mol;
        try {
            mol = searchEngine.smilesToQueryMolecule(query);
        } catch (std::exception &e) {
            LOG(ERROR) << "Cannot parse query " << query << " (" << i + 1 << "): " << e.what();
            exit(1);
        }
        ExperimentInfo info = {
                .stopFlag = false,
                .experimentFinished = false,
                .timeLimit = timeLimit,
        };
        auto checkTimeoutThread = std::thread(checkTimeout, std::ref(info));
        info.startExperiment();
        auto matches = searchEngine.getMatches(*mol, maxResults, info.stopFlag);
        auto experimentDuration = info.finishExperiment();
        info.experimentFinished.test_and_set();
        {
            std::lock_guard<std::mutex> lockGuard(statMutex);
            stat.add(ExperimentStat::Entity{
                    .duration = experimentDuration,
                    .resultsFound = (int) matches.size()
            });
            LOG(INFO) << "Finished " << query << " processing (" << i + 1 << ")"
                      << "\n\tDuration: " << experimentDuration
                      << "\n\tresultsFound: " << matches.size();
        }
        checkTimeoutThread.join();
    });
    statOut << stat << std::endl;
    auto quantiles = stat.quantiles({0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00});
    for (size_t i = 0; i < quantiles.size(); i++) {
        std::cout << (i + 1) * 10 << "%: " << std::setprecision(5) << std::fixed << quantiles[i] << std::endl;
    }
}

int main(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    google::SetLogDestination(google::INFO, "experiment.info");
    FLAGS_alsologtostderr = true;

    ExperimentArgs args(argc, argv);

    auto queries = QueriesParser(args.queriesFile).parse();
    std::ofstream statOut(args.statisticsFile);
    auto smilesDataset = std::make_unique<std::vector<std::string>>();
    std::mutex smilesDatasetMutex;
    SmilesDirParser smilesParser(args.datasetDir);
    smilesParser.parse(*smilesDataset, smilesDatasetMutex);

    if (args.searchEngineType == SearchEngineType::RDKit) {
        auto se = RDKitSearchEngine(std::move(smilesDataset));
        conductExperiment(se, queries, args.maxResults, args.timeLimit, statOut);
    } else if (args.searchEngineType == SearchEngineType::BallTreeRDKit) {
        auto se = BallTreeSearchEngine<RDKitSearchEngine>(std::move(smilesDataset));
        conductExperiment(se, queries, args.maxResults, args.timeLimit, statOut);
    } else if (args.searchEngineType == SearchEngineType::TrieRDKit) {
        LOG(ERROR) << "TrieRDKit is not implemented yet.";
        return 1;
    } else if (args.searchEngineType == SearchEngineType::Indigo) {
        auto se = IndigoSearchEngine(std::move(smilesDataset));
        conductExperiment(se, queries, args.maxResults, args.timeLimit, statOut);
    } else if (args.searchEngineType == SearchEngineType::BallTreeIndigo) {
        auto se = BallTreeSearchEngine<IndigoBruteForceSearchEngine>(std::move(smilesDataset));
        conductExperiment(se, queries, args.maxResults, args.timeLimit, statOut);
    } else if (args.searchEngineType == SearchEngineType::TrieIndigo) {
        auto se = TrieSearchEngine<IndigoBruteForceSearchEngine>(std::move(smilesDataset), 10000);
        conductExperiment(se, queries, args.maxResults, args.timeLimit, statOut);
    } else if (args.searchEngineType == SearchEngineType::IndigoBruteForce) {
        auto se = IndigoBruteForceSearchEngine(std::move(smilesDataset));
        conductExperiment(se, queries, args.maxResults, args.timeLimit, statOut);
    } else {
        LOG(ERROR) << "Specified SearchEngineType is not supported yet";
    }
    return 0;
}
