#pragma once

#include "BallTreeSearchData.h"

namespace qtr {

    class QtrDriveSearchData : public BallTreeSearchData {
    public:

        QtrDriveSearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree,
                           std::shared_ptr<const IdConverter> idConverter, size_t ansCount,
                           size_t threadCount, double timeLimit, bool verificationStage);
    };

} // qtr
