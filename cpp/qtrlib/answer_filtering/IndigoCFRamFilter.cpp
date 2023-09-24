#include "IndigoCFRamFilter.h"

#include <utility>

#include "Utils.h"
#include "molecule/smiles_loader.h"
#include "base_cpp/scanner.h"
#include "IndigoSubstructureMatcher.h"
#include "src/bingo_matcher.h"
#include "Profiling.h"


using namespace qtr;
using namespace indigo;
using namespace indigo_cpp;
using namespace std;
using namespace bingo;


bool IndigoCFRamFilter::operator()(const CIDType &id) {
    ProfileScope("IndigoCFRamFilter");
    try {
        auto molecule = (*_cfStorage)[id];
        molecule->aromatize(AromaticityOptions());
        MoleculeSubstructureMatcher msm(*molecule);
        msm.setQuery(*_queryMolecule);
        return msm.find();
    } catch (const std::exception &e) {
        LOG(ERROR) << e.what();
        return false;
    }
}

IndigoCFRamFilter::IndigoCFRamFilter(std::shared_ptr<CFStorage> cfStorage, const std::string &querySmiles)
        : _cfStorage(std::move(cfStorage)), _querySmiles(querySmiles) {
    BufferScanner scanner(querySmiles.c_str(), (int) querySmiles.size(), false);
    SmilesLoader loader(scanner);
    _queryMolecule = make_unique<QueryMolecule>();
    loader.loadQueryMolecule(*_queryMolecule);
    _queryMolecule->aromatize(AromaticityOptions());
}

unique_ptr<AnswerFilter<CIDType>> IndigoCFRamFilter::copy() {
    return make_unique<IndigoCFRamFilter>(_cfStorage, _querySmiles);
}
