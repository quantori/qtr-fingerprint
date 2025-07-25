#include <glog/logging.h>

#include "ExperimentArgs.h"
#include "io/QueriesParser.h"
#include "io/SmilesDirParser.h"
#include "dataset/SmilesStorage.h"
#include "search/engines/RDKitSearchEngine.h"
#include "search/engines/BallTreeSearchEngine.h"
#include "search/engines/IndigoSearchEngine.h"
#include "benchmarking/BenchmarkRunner.h"
#include "frameworks/IndigoFramework.h"
#include "frameworks/HybridFramework.h"
#include "molecule/query_molecule.h"
#include "utils/Config.h"

int main(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    google::SetLogDestination(google::INFO, "experiment.info");
    FLAGS_alsologtostderr = true;

    try {
        ExperimentArgs args(argc, argv);
        auto querySmiles = QueriesParser(args.queriesFile).parse();
        auto dataSmiles = SmilesDirParser(args.datasetDir).parse();

        Config searchConfig;
        if (args.ballTreeDepth >= 0) {
            searchConfig.set("depth", std::to_string(args.ballTreeDepth));
        }
        searchConfig.set("fpRatio", std::to_string(args.fingerprintRatio));

        BenchmarkArgs benchmarkArgs{
                .queries = querySmiles,
                .maxResults = args.maxResults,
                .timeLimit = args.timeLimit,
                .queriesStatFile = args.queriesStatisticFile,
                .searchEngineStatFile = args.searchEngineStatisticFile,
                .config = searchConfig
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
            case SearchEngineType::BallTreeHybrid: {
                throw std::runtime_error("BallTreeHybrid is not implemented");
//                auto runner = BenchmarkRunner<BallTreeSearchEngine<HybridFramework>>();
//                runner.run(std::move(dataSmiles), benchmarkArgs);
                break;
            }
        }
        ProfilingPool::showStatistics(std::cout);
    } catch (std::exception &e) {
        LOG(ERROR) << e.what();
        return -1;
    } catch (...) {
        LOG(ERROR) << "Something unexpected happened";
        return -1;
    }
    return 0;
}
