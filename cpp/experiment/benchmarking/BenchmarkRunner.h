#pragma once

#include <filesystem>

#include "search/engines/SearchEngineInterface.h"
#include "dataset/SmilesStorage.h"
#include "benchmarking/TimingManager.h"
#include "Profiling.h"
#include "CSVWriter.h"

struct BenchmarkArgs {
    SmilesStorage queries;
    int maxResults;
    double timeLimit;
    std::filesystem::path queriesStatFile;
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
    SearchEngineT buildSearchEngine(SmilesStorage &&dataSmiles) {
        ProfileScope("Total search engine building time");
        return SearchEngineT(std::move(dataSmiles));
    }

    void runQueries(SearchEngineT &searchEngine, const BenchmarkArgs &args) {
        ProfileScope("Total querying running time");
        StatTable queriesStatTable;

        for (size_t i = 0; i < args.queries.size(); i++) {
            auto &smiles = args.queries.smiles(i);
            LOG(INFO) << "Start search [" << i << "]: " << smiles;
            runOneQuery(searchEngine, smiles, args, queriesStatTable);
        }
        // Write statistics
        auto csvTable = queriesStatTable.toCSVTable();
        CSVWriter writer;
        writer.writeToFile(args.queriesStatFile, csvTable);
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

    void runOneQuery(SearchEngineT &searchEngine, const std::string &smiles, const BenchmarkArgs &args,
                     StatTable &statTable) {
        TimingManager timingManager(args.timeLimit);
        timingManager.start();
        SearchQuery query(smiles, args.maxResults, timingManager.getStopFlag());
        try {
            auto result = searchEngine.search(query);
            auto duration = timingManager.finish();
            logOneQueryResult(duration, *result);
            StatRow statRow = result->getStatRow();
            statRow.addEntry("SMILES", smiles);
            statRow.addEntry("FAIL", "No");
            statTable.addRow(statRow);
        } catch (const std::exception &e) {
            LOG(ERROR) << "Query failed: " << e.what();
            StatRow statRow;
            statRow.addEntry("SMILES", smiles);
            statRow.addEntry("FAIL", e.what());
        }

        // TODO: collect per node statistics, write it to the file

        // TODO: collect statistics, write it down
    }
};