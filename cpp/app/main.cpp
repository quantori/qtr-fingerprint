#include "CompleteSearchEngineProfiler.h"
#include "SearchEngineFactory.h"
#include "Utils.h"

#include "IndigoQueryMolecule.h"
#include "indigo.h"

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

#include <glog/logging.h>

using namespace indigo_cpp;
using namespace qtr;

ABSL_FLAG(std::string, path_to_query, "",
          "Path to molecular file to search in database");
ABSL_FLAG(std::string, database_path, "",
          "Path to molecular database");
ABSL_FLAG(std::string, search_engine_type, "bingo",
          "Search engine type: \"bingo\", \"exhaustive\" or \"decision_tree\"");

int main(int argc, char *argv[]) {
    
    google::InitGoogleLogging(argv[0]);
    absl::ParseCommandLine(argc, argv);

    std::string pathToQuery = absl::GetFlag(FLAGS_path_to_query);
    std::string databasePath = absl::GetFlag(FLAGS_database_path);
    std::string searchEngineType = absl::GetFlag(FLAGS_search_engine_type);

    emptyArgument(pathToQuery, "Please specify path_to_query option");
    emptyArgument(databasePath, "Please specify database_path option");
    emptyArgument(searchEngineType, "Please specify search_engine_type option");

    IndigoSessionPtr indigoSessionPtr = IndigoSession::create();
    int queryMoleculeId = indigoLoadQueryMoleculeFromFile(pathToQuery.c_str());
    auto queryMolecule = IndigoQueryMolecule(queryMoleculeId, indigoSessionPtr);

    SearchEngineFactory::SearchEngineType seType = SearchEngineFactory::DECISION_TREE;
    if (searchEngineType == "bingo")
        seType = SearchEngineFactory::BINGO;
    if (searchEngineType == "exhaustive")
        seType = SearchEngineFactory::EXHAUSTIVE;
    std::shared_ptr<SearchEngineInterface> searchEngine = SearchEngineFactory::create(seType, indigoSessionPtr);

    std::vector<IndigoQueryMolecule> queries;
    for (size_t i = 0; i < 10; ++i) {
        queries.emplace_back(queryMolecule);
    }

    CompleteSearchEngineProfiler::profile(databasePath, *searchEngine, queries);
    
    return 0;
}