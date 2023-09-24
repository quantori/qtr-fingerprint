#pragma once

#include <utility>

#include "AnswerFilter.h"
#include "CFStorage.h"
#include "molecule/query_molecule.h"

namespace qtr {

    class IndigoCFRamFilter : public AnswerFilter<CIDType> {
    public:
        IndigoCFRamFilter(std::shared_ptr<CFStorage> cfStorage, const std::string &querySmiles);

        bool operator()(const CIDType &id) override;

        std::unique_ptr<AnswerFilter<CIDType>> copy() override;

    private:
        std::shared_ptr<CFStorage> _cfStorage;
        std::unique_ptr<indigo::QueryMolecule> _queryMolecule;
        std::string _querySmiles;
    };

} // qtr
