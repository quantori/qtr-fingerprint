#include <utility>

#include "IndigoRamFilter.h"
#include "glog/logging.h"
#include "IndigoMolecule.h"
#include "indigo.h"

namespace qtr {

    std::unique_ptr<ByIdAnswerFilter> IndigoRamFilter::copy() {
        return std::make_unique<IndigoRamFilter>(_smilesTable, _querySmiles);
    }

    IndigoRamFilter::IndigoRamFilter(std::shared_ptr<const SmilesTable> smilesTable,
                                     std::shared_ptr<const std::string> querySmiles) :
            IndigoFilter(std::move(querySmiles)), _smilesTable(std::move(smilesTable)) {}

    std::string IndigoRamFilter::getSmiles(CIDType id) {
        return _smilesTable->at(id);
    }

} // namespace qtr
