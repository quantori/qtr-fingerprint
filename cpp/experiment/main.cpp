#include <glog/logging.h>

#include "ExperimentArgs.h"
#include "io/QueriesParser.h"
#include "io/SmilesDirParser.h"
#include "dataset/SmilesStorage.h"
#include "search/engines/RDKit/RDKitSearchEngine.h"
#include "search/engines/BallTree/BallTreeSearchEngine.h"
#include "benchmarking/BenchmarkRunner.h"

int main(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    google::SetLogDestination(google::INFO, "experiment.info");
    FLAGS_alsologtostderr = true;

    ExperimentArgs args(argc, argv);
    auto querySmiles = QueriesParser(args.queriesFile).parse();
    auto dataSmiles = SmilesDirParser(args.datasetDir).parse();
    BenchmarkArgs benchmarkArgs{
            .queries = querySmiles,
            .maxResults = args.maxResults,
            .timeLimit = args.timeLimit,
            .queriesStatFile = args.statisticsFile,
            .searchEngineStatFile = "TODO!"
    };

    switch (args.searchEngineType) {
        case SearchEngineType::BallTreeRDKit: {
            auto runner = BenchmarkRunner<BallTreeSearchEngine<RDKitFramework>>();
            runner.run(std::move(dataSmiles), benchmarkArgs);
            break;
        }
        case SearchEngineType::RDKit: {
            auto runner = BenchmarkRunner<RDKitSearchEngine>();
            runner.run(std::move(dataSmiles), benchmarkArgs);
            break;
        }
        case SearchEngineType::Indigo: {
            throw std::runtime_error("Not implemented error");
            break;
        }
        case SearchEngineType::BallTreeIndigo: {
            throw std::runtime_error("Not implemented error");
            break;
        }
    }
    qtr::ProfilingPool::showStatistics(std::cout);
    return 0;
}
