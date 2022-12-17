#include <utility>

#include "IndigoFilter.h"

#include "glog/logging.h"
#include "IndigoMolecule.h"
#include "indigo.h"

namespace qtr {

    bool IndigoFilter::operator()(CIDType id) {
        const std::string &smiles = _smilesTable->at(id);
        try {
            auto candidateMol = _indigoSessionPtr->loadMolecule(smiles);
            candidateMol.aromatize();
            auto matcher = _indigoSessionPtr->substructureMatcher(candidateMol);
            return bool(indigoMatch(matcher.id(), _queryMolecule.id()));
        }
        catch (std::exception &e) {
            LOG(ERROR) << "Indigo error while filtering. "
                          "Query: " << _querySmiles << ", candidate: " << id << " " << smiles << ", error: "
                       << e.what();
            return false;
        }
    }

    std::unique_ptr<AnswerFilter> IndigoFilter::copy() {
        return std::make_unique<IndigoFilter>(_smilesTable, _querySmiles);
    }

    IndigoFilter::IndigoFilter(std::shared_ptr<const SmilesTable> smilesTable, std::string querySmiles) :
            _smilesTable(std::move(smilesTable)), _querySmiles(std::move(querySmiles)),
            _indigoSessionPtr(indigo_cpp::IndigoSession::create()),
            _queryMolecule(_indigoSessionPtr->loadQueryMolecule(_querySmiles)) {
        _queryMolecule.aromatize();
    }

} // namespace qtr
