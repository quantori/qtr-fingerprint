#pragma once

#include <memory>
#include <vector>

#include "CFStorage.h"

#include "SmilesTable.h"
#include "QtrSearchData.h"

namespace qtr {

    class QtrRamSearchData : public QtrSearchData {
    public:

        QtrRamSearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree,
                         std::shared_ptr<const IdConverter> idConverter, size_t ansCount,
                         size_t threadCount, double timeLimit,
                         std::shared_ptr<CFStorage> cfStorage,
                         std::shared_ptr<const std::vector<PropertiesFilter::Properties>> propertiesTable);


        std::shared_ptr<CFStorage> cfStorage;
        std::shared_ptr<const std::vector<PropertiesFilter::Properties>> propertiesTable;
    };

} // qtr
