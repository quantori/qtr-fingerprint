#pragma once

#include <string>
#include <atomic>

#include "IndigoFilter.h"
#include "SmilesTable.h"

#include "IndigoQueryMolecule.h"
#include "IndigoSubstructureMatcher.h"

namespace qtr {

    class IndigoRamSmilesFilter : public IndigoFilter {
    public:
        IndigoRamSmilesFilter(std::shared_ptr<const SmilesTable> smilesTable, std::shared_ptr<const std::string> querySmiles);

        std::unique_ptr<AnswerFilter> copy() override;

    private:
        std::shared_ptr<indigo_cpp::IndigoMolecule> getMolecule(CIDType id) override;

        std::shared_ptr<const SmilesTable> _smilesTable;
    };

} // namespace qtr
