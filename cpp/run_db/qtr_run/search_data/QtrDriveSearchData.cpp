#include "QtrDriveSearchData.h"

#include <utility>

namespace qtr {
    QtrDriveSearchData::QtrDriveSearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree,
                                           std::shared_ptr<const IdConverter> idConverter, TimeTicker &timeTicker,
                                           size_t ansCount, size_t threadCount, double timeLimit) :
            QtrSearchData(std::move(ballTree), std::move(idConverter), timeTicker, ansCount, threadCount, timeLimit) {}
} // qtr