#include <memory>
#include <utility>

#include "IndigoRamMoleculeFilter.h"
#include "IndigoMolecule.h"


using namespace std;

namespace qtr {
    IndigoRamMoleculeFilter::IndigoRamMoleculeFilter(
            shared_ptr<const map <CIDType, indigo_cpp::IndigoMolecule>> moleculesTable,
            shared_ptr<const string> querySmiles) :
            IndigoFilter(std::move(querySmiles)),
            _moleculesTable(std::move(moleculesTable)) {}

    unique_ptr <AnswerFilter> IndigoRamMoleculeFilter::copy() {
        return make_unique<IndigoRamMoleculeFilter>(_moleculesTable, _querySmiles);
    }

    shared_ptr <indigo_cpp::IndigoMolecule> IndigoRamMoleculeFilter::getMolecule(CIDType id) {
        return make_shared<indigo_cpp::IndigoMolecule>(_moleculesTable->at(id));
    }

} // qtr