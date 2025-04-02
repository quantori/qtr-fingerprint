#pragma once

#include <filesystem>

#include "search/engines/SearchEngineInterface.h"
#include "dataset/SmilesStorage.h"
#include "benchmarking/TimingManager.h"
#include "Profiling.h"
#include "benchmarking/stats/QueryStats.h"

struct BenchmarkArgs {
    SmilesStorage queries;
    int maxResults;
    double timeLimit;
    std::filesystem::path queriesStatFile; // TODO: it is not used now
    std::filesystem::path searchEngineStatFile; // TODO: it is not used now
};

template<typename SearchEngineT> requires SearchEngineInterface<SearchEngineT>
class BenchmarkRunner {
public:
    using FrameworkT = SearchEngineT::FrameworkT;

    BenchmarkRunner() = default;

    ~BenchmarkRunner() {
        // TODO: print quick summary on exit
    }

    void run(SmilesStorage &&dataSmiles, BenchmarkArgs &args) {
        ProfileScope("Total benchmark running time");
        auto searchEngine = buildSearchEngine(std::move(dataSmiles));
        runQueries(searchEngine, args);
    }

private:
//    QueryStats _queryStats; // todo

    SearchEngineT buildSearchEngine(SmilesStorage &&dataSmiles) {
        ProfileScope("Total search engine building time");
        return SearchEngineT(std::move(dataSmiles));
    }

    void runQueries(SearchEngineT &searchEngine, const BenchmarkArgs &args) {
        ProfileScope("Total querying running time");
        for (size_t i = 0; i < args.queries.size(); i++) {
            auto &smiles = args.queries.smiles(i);
            LOG(INFO) << "Start search [" << i << "]: " << smiles;
            runOneQuery(searchEngine, smiles, args);
        }
    }

    void runOneQuery(SearchEngineT &searchEngine, const std::string &smiles, const BenchmarkArgs &args) {
        TimingManager timingManager(args.timeLimit);
        timingManager.start();
        SearchQuery query(smiles, args.maxResults, timingManager.getStopFlag());
        auto result = searchEngine.search(query);
        auto duration = timingManager.finish();
        LOG(INFO) << "Search finished in " << duration << " seconds. Found " << result->size() << " answers";
        // TODO: collect statistics, write it down
    }
};