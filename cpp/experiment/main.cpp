#include <glog/logging.h>

#include "ExperimentArgs.h"
#include "io/QueriesParser.h"
#include "io/SmilesDirParser.h"
#include "dataset/SmilesStorage.h"
#include "search/engines/RDKit/RDKitSearchEngine.h"
#include "search/engines/BallTree/BallTreeSearchEngine.h"
#include "search/engines/Indigo/IndigoSearchEngine.h"
#include "benchmarking/BenchmarkRunner.h"

int main(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    google::SetLogDestination(google::INFO, "experiment.info");
    FLAGS_alsologtostderr = true;

    try {
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
                auto runner = BenchmarkRunner<IndigoSearchEngine>();
                runner.run(std::move(dataSmiles), benchmarkArgs);
                break;
            }
            case SearchEngineType::BallTreeIndigo: {
                auto runner = BenchmarkRunner<BallTreeSearchEngine<IndigoFramework>>();
                runner.run(std::move(dataSmiles), benchmarkArgs);
                break;
            }
        }
        qtr::ProfilingPool::showStatistics(std::cout);
    } catch (std::exception &e) {
        LOG(ERROR) << e.what();
        return -1;
    } catch (...) {
        LOG(ERROR) << "Something unexpected happened";
        return -1;
    }
    return 0;
}
