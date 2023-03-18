#include <utility>
#include <memory>

#include "IndigoRamSmilesFilter.h"
#include "glog/logging.h"
#include "IndigoMolecule.h"
#include "indigo.h"

using namespace std;

namespace qtr {
    unique_ptr<AnswerFilter> IndigoRamSmilesFilter::copy() {
        return make_unique<IndigoRamSmilesFilter>(_smilesTable, _querySmiles);
    }

    IndigoRamSmilesFilter::IndigoRamSmilesFilter(shared_ptr<const SmilesTable> smilesTable,
                                     shared_ptr<const string> querySmiles) :
            IndigoFilter(std::move(querySmiles)), _smilesTable(std::move(smilesTable)) {}

    shared_ptr<indigo_cpp::IndigoMolecule> IndigoRamSmilesFilter::getMolecule(CIDType id) {
        auto mol = std::make_shared<indigo_cpp::IndigoMolecule>(_indigoSessionPtr->loadMolecule(_smilesTable->at(id)));
        mol->aromatize();
        return mol;
    }

} // namespace qtr
