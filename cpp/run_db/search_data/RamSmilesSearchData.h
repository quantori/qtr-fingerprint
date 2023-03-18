#pragma once

#include <memory>
#include <vector>

#include "SearchData.h"
#include "SmilesTable.h"

namespace qtr {

    class RamSmilesSearchData : public SearchData {
    public:

        RamSmilesSearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree,
                            std::shared_ptr<const IdConverter> idConverter, TimeTicker &timeTicker, size_t ansCount,
                            size_t threadCount, std::shared_ptr<const SmilesTable> smilesTable,
                            std::shared_ptr<const std::vector<PropertiesFilter::Properties>> propertiesTable);

        [[nodiscard]] DbType getDbType() const override;

        std::shared_ptr<const SmilesTable> smilesTable;
        std::shared_ptr<const std::vector<PropertiesFilter::Properties>> propertiesTable;
    };

} // qtr
