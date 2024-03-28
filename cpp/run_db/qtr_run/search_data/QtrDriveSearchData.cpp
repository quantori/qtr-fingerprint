#include "QtrDriveSearchData.h"

#include <utility>

namespace qtr {
    QtrDriveSearchData::QtrDriveSearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree,
                                           std::shared_ptr<const IdConverter> idConverter, size_t ansCount,
                                           size_t threadCount, double timeLimit, bool verificationStage) :
            BallTreeSearchData(std::move(ballTree), std::move(idConverter), ansCount, threadCount, timeLimit,
                               verificationStage) {}
} // qtr