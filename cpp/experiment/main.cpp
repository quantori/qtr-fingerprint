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

#include "Profiling.h"

struct ExperimentInfo {
    bool stopFlag;
    std::atomic<bool> experimentFinished{false};
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
            LOG(INFO) << "Timeout detected";
        }
    }

    void startExperiment() {
        std::lock_guard<std::mutex> lockGuard(mutex);
        stopFlag = false;
        startTime = std::chrono::high_resolution_clock::now();
    }

    [[nodiscard]] double finishExperiment() {
        std::chrono::duration<double> duration = std::chrono::high_resolution_clock::now() - startTime;
        experimentFinished.store(true, std::memory_order_release);
        return duration.count();
    }
};

void checkTimeout(ExperimentInfo &info) {
    std::chrono::duration<double> sleepDuration(std::max(0.001, std::min(0.01, info.timeLimit)));
    while (!info.experimentFinished.load(std::memory_order_relaxed)) {
        std::this_thread::sleep_for(sleepDuration);
        info.checkTimeout();
        if (info.stopFlag) {
            LOG(INFO) << "Timeout thread detected stop flag, waiting for experiment to finish...";
        }
    }
    LOG(INFO) << "Timeout thread finishing...";
}

class ExperimentSession {
    ExperimentInfo info;
    std::jthread timeoutThread;

public:
    explicit ExperimentSession(double timeLimit) 
        : info{
            .stopFlag = false,
            .experimentFinished = false,
            .timeLimit = timeLimit
          }
        , timeoutThread(checkTimeout, std::ref(info))
    {
    }

    ~ExperimentSession() {
        info.experimentFinished.store(true, std::memory_order_release);
    }

    void start() {
        info.startExperiment();
    }

    [[nodiscard]] double finish() {
        return info.finishExperiment();
    }

    [[nodiscard]] bool& stopFlag() {
        return info.stopFlag;
    }
};

template<SearchEngine SE>
void conductExperiment(SE &searchEngine,
                       const std::vector<std::string> &queries,
                       int maxResults,
                       double timeLimit,
                       std::ofstream &statOut) {
    ProfileScope("Experiment");
    ExperimentStat stat;
    std::mutex statMutex;
    std::atomic_int64_t counter = 0;
    
    std::for_each(std::execution::seq, queries.begin(), queries.end(), [&](const std::string &query) {
        auto i = counter++;
        LOG(INFO) << "Start " << query << " processing (" << i + 1 << ")";
        
        try {
            decltype(searchEngine.smilesToQueryMolecule(query)) mol;
            try {
                mol = searchEngine.smilesToQueryMolecule(query);
            } catch (std::exception &e) {
                LOG(ERROR) << "Cannot parse query " << query << " (" << i + 1 << "): " << e.what();
                return;
            }

            ExperimentSession session(timeLimit);
            session.start();
            
            LOG(INFO) << "Starting getMatches for query " << i + 1;
            auto matches = searchEngine.getMatches(*mol, maxResults, session.stopFlag());
            LOG(INFO) << "getMatches completed for query " << i + 1;
            
            auto experimentDuration = session.finish();
            
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
        } catch (const std::exception& e) {
            LOG(ERROR) << "Exception in experiment loop: " << e.what() 
                      << " for query " << query << " (" << i + 1 << ")";
        } catch (...) {
            LOG(ERROR) << "Unknown exception in experiment loop"
                      << " for query " << query << " (" << i + 1 << ")";
        }
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
        LOG(ERROR) << "TrieIndigo is not implemented yet.";
        return 1;
//        auto se = TrieSearchEngine<IndigoBruteForceSearchEngine>(std::move(smilesDataset), 10000);
//        conductExperiment(se, queries, args.maxResults, args.timeLimit, statOut);
    } else if (args.searchEngineType == SearchEngineType::IndigoBruteForce) {
        auto se = IndigoBruteForceSearchEngine(std::move(smilesDataset));
        conductExperiment(se, queries, args.maxResults, args.timeLimit, statOut);
    } else {
        LOG(ERROR) << "Specified SearchEngineType is not supported yet";
    }
    qtr::ProfilingPool::showStatistics(std::cout);
    return 0;
}
