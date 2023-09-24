#pragma once

#include <atomic>
#include <string>

#include "AnswerFilter.h"
#include "SmilesTable.h"
#include "BallTreeTypes.h"

#include "IndigoQueryMolecule.h"
#include "IndigoSubstructureMatcher.h"

namespace qtr {

    class IndigoSmilesFilter : public ByIdAnswerFilter {
    public:
        bool operator()(const CIDType &id) override;

        explicit IndigoSmilesFilter(std::shared_ptr<const std::string> querySmiles);

    protected:
        virtual std::string getSmiles(CIDType id) = 0;

        std::shared_ptr<const std::string> _querySmiles;
        indigo_cpp::IndigoSessionPtr _indigoSessionPtr;
        indigo_cpp::IndigoQueryMolecule _queryMolecule;
    };

} // qtr
