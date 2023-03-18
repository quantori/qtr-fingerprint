#pragma once

#include <memory>
#include <map>

#include "IndigoFilter.h"

namespace qtr {

    class IndigoRamMoleculeFilter : public IndigoFilter {
    public:
        IndigoRamMoleculeFilter(std::shared_ptr<const std::map<CIDType, indigo_cpp::IndigoMolecule>> moleculesTable,
                                std::shared_ptr<const std::string> querySmiles);

        std::unique_ptr<AnswerFilter> copy() override;

    private:
        std::shared_ptr<indigo_cpp::IndigoMolecule> getMolecule(CIDType id) override;

        std::shared_ptr<const std::map<CIDType, indigo_cpp::IndigoMolecule>> _moleculesTable;
    };

} // qtr

