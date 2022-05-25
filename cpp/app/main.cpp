#include "CompleteSearchEngineProfiler.h"
#include "SearchEngineFactory.h"
#include "IndigoQueryMolecule.h"
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

#include <glog/logging.h>

using namespace indigo_cpp;

ABSL_FLAG(std::string, path_to_query, "",
          "Path to molecular file to search in database");
ABSL_FLAG(std::string, database_path, "",
          "Path to molecular database");

int main(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    absl::ParseCommandLine(argc, argv);
    std::string pathToQuery = absl::GetFlag(FLAGS_path_to_query);
    std::string databasePath = absl::GetFlag(FLAGS_database_path);
    emptyArgument(pathToQuery, "Please specify path_to_query option");
    emptyArgument(databasePath, "Please specify database_path option");
    IndigoSessionPtr indigoSessionPtr = IndigoSession::create();
    int queryMoleculeId = indigoLoadQueryMoleculeFromFile(pathToQuery.c_str());
    auto queryMolecule = IndigoQueryMolecule(queryMoleculeId, indigoSessionPtr);
    std::shared_ptr<SearchEngineInterface> searchEngine = SearchEngineFactory::create(indigoSessionPtr);
    std::vector<IndigoQueryMolecule> queries;
    for (size_t i = 0; i < 10; ++i) {
        queries.emplace_back(queryMolecule);
    }
    SearchEngineProfiler p(databasePath, *searchEngine);
    CompleteSearchEngineProfiler::profile(databasePath, *searchEngine, queries);
    return 0;
}