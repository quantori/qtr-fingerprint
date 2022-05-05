#include "CompleteSearchEngineProfiler.h"
#include "MockSearchEngine.h"
#include "IndigoQueryMolecule.h"

using namespace indigo_cpp;

int main(void) {
    IndigoSessionPtr indigoSessionPtr = IndigoSession::create();
    IndigoMolecule molecule = indigoSessionPtr->loadMoleculeFromFile("query.mol");
    IndigoQueryMolecule queryMolecule = IndigoQueryMolecule(molecule.id(), indigoSessionPtr);
    std::string pathToFile = "1.sdf.gz";
    MockSearchEngine searchEngine;
    std::vector<IndigoQueryMolecule> queries;
    for (size_t i = 0; i < 100; ++i) {
        queries.emplace_back(queryMolecule);
    }
    SearchEngineProfiler p(pathToFile, searchEngine);
    CompleteSearchEngineProfiler::profile(pathToFile, searchEngine, queries);
    return 0;
}