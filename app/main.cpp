#include "CompleteSearchEngineProfiler.h"
#include "SearchEngineFactory.h"
#include "IndigoQueryMolecule.h"

using namespace indigo_cpp;

qword sessionId;

int main(void) {
    IndigoSessionPtr indigoSessionPtr = IndigoSession::create();
    sessionId = indigoSessionPtr->getSessionId();
    int queryMoleculeId = indigoLoadQueryMoleculeFromFile("query.mol");
    auto queryMolecule = IndigoQueryMolecule(queryMoleculeId, indigoSessionPtr);
    std::string pathToFile = "119697";
    std::shared_ptr<SearchEngineInterface> searchEngine = SearchEngineFactory::create();
    std::vector<IndigoQueryMolecule> queries;
    for (size_t i = 0; i < 1; ++i) {
        queries.emplace_back(queryMolecule);
    }
    SearchEngineProfiler p(pathToFile, *searchEngine);
    CompleteSearchEngineProfiler::profile(pathToFile, *searchEngine, queries);
    return 0;
}