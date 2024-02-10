#include "QtrRamSearchData.h"

#include <utility>

namespace qtr {
    QtrRamSearchData::QtrRamSearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree,
                                       std::shared_ptr<const IdConverter> idConverter, size_t ansCount,
                                       size_t threadCount, double timeLimit, std::shared_ptr<CFStorage> cfStorage,
                                       std::shared_ptr<const std::vector<PropertiesFilter::Properties>> propertiesTable,
                                       bool verificationStage)
            : QtrSearchData(std::move(ballTree), std::move(idConverter), ansCount, threadCount, timeLimit,
                            verificationStage), cfStorage(std::move(cfStorage)),
              propertiesTable(std::move(propertiesTable)) {
    }
} // qtr