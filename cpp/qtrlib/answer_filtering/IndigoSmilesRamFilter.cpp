#include <utility>

#include "IndigoSmilesRamFilter.h"
#include "glog/logging.h"
#include "IndigoMolecule.h"
#include "indigo.h"

namespace qtr {

    std::unique_ptr<ByIdAnswerFilter> IndigoSmilesRamFilter::copy() {
        return std::make_unique<IndigoSmilesRamFilter>(_smilesTable, _querySmiles);
    }

    IndigoSmilesRamFilter::IndigoSmilesRamFilter(std::shared_ptr<const SmilesTable> smilesTable,
                                                 std::shared_ptr<const std::string> querySmiles) :
            IndigoSmilesFilter(std::move(querySmiles)), _smilesTable(std::move(smilesTable)) {}

    std::string IndigoSmilesRamFilter::getSmiles(CIDType id) {
        return _smilesTable->at(id);
    }

} // namespace qtr
