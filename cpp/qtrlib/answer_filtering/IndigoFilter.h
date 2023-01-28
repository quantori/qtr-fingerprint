#pragma once

#include <string>

#include "AnswerFilter.h"
#include "SmilesTable.h"

#include "IndigoQueryMolecule.h"
#include "IndigoSubstructureMatcher.h"

namespace qtr {

    class IndigoFilter : public AnswerFilter {
    public:
        IndigoFilter(std::shared_ptr<const SmilesTable> smilesTable, std::shared_ptr<const std::string> querySmiles);

        bool operator()(CIDType id) override;

        std::unique_ptr<AnswerFilter> copy() override;

    private:
        std::shared_ptr<const SmilesTable> _smilesTable;
        std::shared_ptr<const std::string> _querySmiles;
        indigo_cpp::IndigoSessionPtr _indigoSessionPtr;
        indigo_cpp::IndigoQueryMolecule _queryMolecule;
    };

} // namespace qtr