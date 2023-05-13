#pragma once

#include <memory>
#include <vector>

#include "SearchData.h"
#include "SmilesTable.h"

namespace qtr {

    class RamSearchData : public SearchData {
    public:

        RamSearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree,
                      std::shared_ptr<const IdConverter> idConverter, TimeTicker &timeTicker, size_t ansCount,
                      size_t threadCount, double timeLimit, std::shared_ptr<const SmilesTable> smilesTable,
                      std::shared_ptr<const std::vector<PropertiesFilter::Properties>> propertiesTable);

        [[nodiscard]] DerivedClasses getClass() const override;


        std::shared_ptr<const SmilesTable> smilesTable;
        std::shared_ptr<const std::vector<PropertiesFilter::Properties>> propertiesTable;
    };

} // qtr
