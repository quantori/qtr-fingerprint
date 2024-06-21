#pragma once

#include "answer_filtering/AnswerFilter.h"
#include "GraphMol/SubstructLibrary/SubstructLibrary.h"

namespace qtr {

    class RDKitSubstructFilter : public AnswerFilter<CIDType> {
    public:
        RDKitSubstructFilter(std::shared_ptr<RDKit::MolHolderBase> molHolder, const std::string &querySmiles);

        RDKitSubstructFilter(std::shared_ptr<RDKit::MolHolderBase> molHolder,
                             std::shared_ptr<RDKit::ROMol> queryMolecule);

        std::unique_ptr<AnswerFilter<CIDType>> copy() override;

        bool operator()(const CIDType& id) override;

    private:
        std::shared_ptr<RDKit::MolHolderBase> _molHolder;
        std::shared_ptr<RDKit::ROMol> _queryMolecule;
    };

} // qtr

