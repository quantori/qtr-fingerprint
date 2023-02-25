#include "RamSearchData.h"

#include <utility>

namespace qtr {
    SearchData::DerivedClasses RamSearchData::getClass() const {
        return SearchData::DerivedClasses::RamSearchData;
    }

    RamSearchData::RamSearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree,
                                 std::shared_ptr<const IdConverter> idConverter, TimeTicker &timeTicker,
                                 size_t ansCount,
                                 size_t threadCount, std::shared_ptr<const SmilesTable> smilesTable,
                                 std::shared_ptr<const std::vector<PropertiesFilter::Properties>> propertiesTable) :
            SearchData(std::move(ballTree), std::move(idConverter), timeTicker, ansCount, threadCount),
            smilesTable(std::move(smilesTable)), propertiesTable(std::move(propertiesTable)) {}
} // qtr