#include "IndigoCFRamFilter.h"

#include <utility>

#include "Utils.h"
#include "molecule/smiles_loader.h"
#include "base_cpp/scanner.h"
#include "IndigoSubstructureMatcher.h"
#include "src/bingo_matcher.h"


using namespace qtr;
using namespace indigo;
using namespace indigo_cpp;
using namespace std;
using namespace bingo;


bool IndigoCFRamFilter::operator()(const CIDType &id) {
    try {
        auto molecule = (*_cfStorage)[id];
        MoleculeSubstructureMatcher msm(*molecule);
        msm.setQuery(*_queryMolecule);
        return msm.find();

    } catch (const std::exception &e) {
        logErrorAndExit(e.what());
    }
}

IndigoCFRamFilter::IndigoCFRamFilter(std::shared_ptr<CFStorage> cfStorage, const std::string &querySmiles)
        : _cfStorage(std::move(cfStorage)) {
    BufferScanner scanner(querySmiles.c_str(), querySmiles.size(), false);
    SmilesLoader loader(scanner);
    loader.loadQueryMolecule(*_queryMolecule);
}

unique_ptr<AnswerFilter<CIDType>> IndigoCFRamFilter::copy() {
    return make_unique<IndigoCFRamFilter>(*this);
}

IndigoCFRamFilter::IndigoCFRamFilter(const IndigoCFRamFilter &other) : _cfStorage(other._cfStorage),
                                                                       _queryMolecule(other._queryMolecule) {

}
// qtr