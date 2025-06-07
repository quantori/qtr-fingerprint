#pragma once

#include <filesystem>

#include "search/engines/SearchEngineInterface.h"
#include "dataset/SmilesStorage.h"
#include "benchmarking/TimingManager.h"
#include "Profiling.h"
#include "CSVWriter.h"
#include "utils/Config.h"

struct BenchmarkArgs {
    SmilesStorage queries;
    int maxResults;
    double timeLimit;
    std::filesystem::path queriesStatFile;
    std::filesystem::path searchEngineStatFile;
    Config config;
};

template<typename SearchEngineT> requires SearchEngineInterface<SearchEngineT>
class BenchmarkRunner {
public:
    using FrameworkT = SearchEngineT::FrameworkT;

    BenchmarkRunner() = default;

    void run(SmilesStorage &&dataSmiles, BenchmarkArgs &args) {
        ProfileScope("Total benchmark running time");
        auto framework = FrameworkT(args.config);
        auto searchEngine = buildSearchEngine(std::move(framework), std::move(dataSmiles), args.config);
        runQueries(searchEngine, args);
    }

private:
    SearchEngineT
    buildSearchEngine(FrameworkT framework, SmilesStorage &&dataSmiles, const Config &config) {
        ProfileScope("Total search engine building time");
        return SearchEngineT(std::move(framework), std::move(dataSmiles), config);
    }

    void runQueries(SearchEngineT &searchEngine, const BenchmarkArgs &args) {
        ProfileScope("Total querying running time");
        StatTable queriesStatTable;

        for (size_t i = 0; i < args.queries.size(); i++) {
            auto &smiles = args.queries.smiles(i);
            LOG(INFO) << "Start search [" << i << "]: " << smiles;
            runOneQuery(searchEngine, smiles, args, queriesStatTable);
        }
        writeStatistics(queriesStatTable, searchEngine, args);
    }

    void logOneQueryResult(double duration, const SearchResult<typename SearchEngineT::ResultT> &result) {
        LOG(INFO) << "Search finished in " << duration << " seconds. Found " << result.size() << " answers";
        if (result.size() > 0) {
            LOG(INFO) << "Example results:";
            for (size_t i = 0; i < result.size() && i < 5; i++) {
                const auto &res = result.get(i);
                if constexpr (std::is_same_v<typename FrameworkT::MoleculeT, typename SearchEngineT::ResultT>) {
                    auto resSmiles = FrameworkT::moleculeToSmiles(res);
                    LOG(INFO) << resSmiles;
                } else {
                    LOG(INFO) << res;
                }
            }
        }
    }

    void
    writeStatistics(const StatTable &queriesStatTable, const SearchEngineT &searchEngine, const BenchmarkArgs &args) {
        CSVWriter writer;

        auto queriesCSVTable = queriesStatTable.toCSVTable();
        writer.writeToFile(args.queriesStatFile, queriesCSVTable);


        StatTable searchEngineStatTable = searchEngine.getStat();
        auto searchEngineCSVTable = searchEngineStatTable.toCSVTable();
        writer.writeToFile(args.searchEngineStatFile, searchEngineCSVTable);
    }

    void runOneQuery(SearchEngineT &searchEngine, const std::string &smiles, const BenchmarkArgs &args,
                     StatTable &statTable) {
        TimingManager timingManager(args.timeLimit);
        timingManager.start();
        SearchQuery query(smiles, args.maxResults, timingManager.getStopFlag());
        try {
            std::unique_ptr<SearchResult<typename SearchEngineT::ResultT>> result = searchEngine.search(query);
            auto duration = timingManager.finish();
            logOneQueryResult(duration, *result);
            StatRow statRow = result->getStatRow();
            statRow.addEntry("SMILES", smiles);
            statRow.addEntry("FAIL", "No");
            statRow.addEntry("ANS_COUNT", result->size());
            statRow.addEntry("DURATION", duration);
            statTable.addRow(statRow);
        } catch (const std::exception &e) {
            LOG(ERROR) << "Query failed: " << e.what();
            StatRow statRow;
            statRow.addEntry("SMILES", smiles);
            statRow.addEntry("FAIL", e.what());
            statTable.addRow(statRow);
        }
    }
};