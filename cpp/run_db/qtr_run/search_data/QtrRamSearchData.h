#pragma once

#include <memory>
#include <vector>

#include "CFStorage.h"

#include "SmilesTable.h"
#include "BallTreeSearchData.h"

namespace qtr {

    class QtrRamSearchData : public BallTreeSearchData {
    public:

        QtrRamSearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree,
                         std::shared_ptr<const IdConverter> idConverter, size_t ansCount,
                         size_t threadCount, double timeLimit, std::shared_ptr<CFStorage> cfStorage,
                         std::shared_ptr<const std::vector<PropertiesFilter::Properties>> propertiesTable,
                         bool verificationStage);


        std::shared_ptr<CFStorage> cfStorage;
        std::shared_ptr<const std::vector<PropertiesFilter::Properties>> propertiesTable;
    };

} // qtr
