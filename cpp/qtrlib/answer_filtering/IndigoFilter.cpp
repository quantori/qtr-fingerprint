#include <utility>

#include "IndigoFilter.h"

#include "glog/logging.h"
#include "IndigoMolecule.h"
#include "indigo.h"

namespace qtr {

    bool IndigoFilter::operator()(CIDType id) {
        auto startTime = std::chrono::high_resolution_clock::now();
        const std::string &smiles = _smilesTable->at(id);
        bool result;
        try {
            auto candidateMol = _indigoSessionPtr->loadMolecule(smiles);
            candidateMol.aromatize();
            auto matcher = _indigoSessionPtr->substructureMatcher(candidateMol);
            result = bool(indigoMatch(matcher.id(), _queryMolecule.id()));
        }
        catch (std::exception &e) {
            LOG(ERROR) << "Indigo error while filtering. "
                          "Query: " << _querySmiles << ", candidate: " << id << " " << smiles << ", error: "
                       << e.what();
            result = false;
        }
        std::chrono::duration<double> duration = std::chrono::high_resolution_clock::now() - startTime;
        _timer += duration.count();
        return result;
    }

    std::unique_ptr<AnswerFilter> IndigoFilter::copy() {
        return std::make_unique<IndigoFilter>(_smilesTable, _querySmiles);
    }

    IndigoFilter::IndigoFilter(std::shared_ptr<const SmilesTable> smilesTable,
                               std::shared_ptr<const std::string> querySmiles) :
            _smilesTable(std::move(smilesTable)), _querySmiles(std::move(querySmiles)),
            _indigoSessionPtr(indigo_cpp::IndigoSession::create()),
            _queryMolecule(_indigoSessionPtr->loadQueryMolecule(*_querySmiles)) {
        _queryMolecule.aromatize();
    }

    IndigoFilter::~IndigoFilter() {
        indigoFilteringTimer += _timer;
    }

} // namespace qtr
