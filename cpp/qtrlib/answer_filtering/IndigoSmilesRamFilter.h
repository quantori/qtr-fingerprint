#pragma once

#include <string>
#include <atomic>

#include "IndigoSmilesFilter.h"
#include "SmilesTable.h"

#include "IndigoQueryMolecule.h"
#include "IndigoSubstructureMatcher.h"

namespace qtr {


    class IndigoSmilesRamFilter : public IndigoSmilesFilter {
    public:
        IndigoSmilesRamFilter(std::shared_ptr<const SmilesTable> smilesTable, std::shared_ptr<const std::string> querySmiles);

        std::unique_ptr<AnswerFilter> copy() override;

    private:
        std::string getSmiles(CIDType id) override;

        std::shared_ptr<const SmilesTable> _smilesTable;
    };

} // namespace qtr
