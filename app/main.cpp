#include "CompleteSearchEngineProfiler.h"
#include "MockSearchEngine.h"

using namespace indigo_cpp;

int main(void) {
    IndigoSessionPtr indigoSessionPtr = IndigoSession::create();
    auto queryMolecule = indigoSessionPtr->loadMoleculeFromFile("query.mol");
    std::string pathToFile = "1.sdf.gz";
    MockSearchEngine searchEngine;
    std::vector<IndigoMolecule> queries;
    for (size_t i = 0; i < 100; ++i) {
        queries.emplace_back(queryMolecule);
    }
    SearchEngineProfiler p(pathToFile, searchEngine);
    CompleteSearchEngineProfiler::profile(pathToFile, searchEngine,queries);
    return 0;
}