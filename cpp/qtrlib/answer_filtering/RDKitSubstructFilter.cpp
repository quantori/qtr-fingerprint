#include "RDKitSubstructFilter.h"

using namespace std;

namespace qtr {
    RDKitSubstructFilter::RDKitSubstructFilter(std::shared_ptr<RDKit::MolHolderBase> molHolder,
                                               const std::string &querySmiles) : _molHolder(std::move(molHolder)) {
        try {
            _queryMolecule = shared_ptr<RDKit::ROMol>(RDKit::SmilesToMol(querySmiles));
            if (_queryMolecule == nullptr) {
                throw runtime_error("Cannot create RDKitSubstructFilter for " + querySmiles);
            }
        } catch (const exception &e) {
            throw runtime_error("Cannot create RDKitSubstructFilter for " + querySmiles + " because of " + e.what());
        }
    }

    RDKitSubstructFilter::RDKitSubstructFilter(std::shared_ptr<RDKit::MolHolderBase> molHolder,
                                               std::shared_ptr<RDKit::ROMol> queryMolecule) : _molHolder(
            std::move(molHolder)), _queryMolecule(std::move(queryMolecule)) {
    }

    unique_ptr <AnswerFilter<CIDType>> RDKitSubstructFilter::copy() {
        return make_unique<RDKitSubstructFilter>(_molHolder, _queryMolecule);
    }

    bool RDKitSubstructFilter::operator()(const CIDType &id) {
        auto mol = _molHolder->getMol(id);
        RDKit::SubstructMatchParameters params;
        params.recursionPossible = true;
        params.useChirality = true;
        params.useQueryQueryMatches = false;
        auto matches = RDKit::SubstructMatch(*mol, *_queryMolecule, params);
        return !matches.empty();
    }
} // qtr