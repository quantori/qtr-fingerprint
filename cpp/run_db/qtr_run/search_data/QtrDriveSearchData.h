#pragma once

#include "QtrSearchData.h"

namespace qtr {

    class QtrDriveSearchData : public QtrSearchData {
    public:

        QtrDriveSearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree,
                           std::shared_ptr<const IdConverter> idConverter,
                           TimeTicker &timeTicker, size_t ansCount,
                           size_t threadCount, double timeLimit);
    };

} // qtr
