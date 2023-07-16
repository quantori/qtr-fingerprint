#include "QtrRamSearchData.h"

#include <utility>

namespace qtr {
    QtrRamSearchData::QtrRamSearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree,
                                       std::shared_ptr<const IdConverter> idConverter, TimeTicker &timeTicker,
                                       size_t ansCount, size_t threadCount, double timeLimit,
                                       std::shared_ptr<const SmilesTable> smilesTable,
                                       std::shared_ptr<const std::vector<PropertiesFilter::Properties>> propertiesTable)
            : QtrSearchData(std::move(ballTree), std::move(idConverter), timeTicker, ansCount, threadCount, timeLimit),
              smilesTable(std::move(smilesTable)), propertiesTable(std::move(propertiesTable)) {}
} // qtr