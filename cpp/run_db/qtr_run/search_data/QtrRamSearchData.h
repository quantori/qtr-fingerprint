#pragma once

#include <memory>
#include <vector>

#include "CFStorage.h"

#include "SmilesTable.h"
#include "BallTreeSearchData.h"

#include <GraphMol/SubstructLibrary/SubstructLibrary.h>

namespace qtr {

    class QtrRamSearchData : public BallTreeSearchData {
    public:

        QtrRamSearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree,
                         std::shared_ptr<const IdConverter> idConverter, size_t ansCount,
                         size_t threadCount, double timeLimit, std::shared_ptr<CFStorage> cfStorage,
                         std::shared_ptr<RDKit::MolHolderBase> molHolder,
                         std::shared_ptr<const std::vector<PropertiesFilter::Properties>> propertiesTable,
                         bool verificationStage);

        BaseLibrary getBaseLibrary() const override;

        std::shared_ptr<CFStorage> cfStorage;
        std::shared_ptr<RDKit::MolHolderBase> molHolder;
        std::shared_ptr<const std::vector<PropertiesFilter::Properties>> propertiesTable;
    };

} // qtr
