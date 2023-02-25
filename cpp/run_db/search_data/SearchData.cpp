#include "SearchData.h"

#include <utility>

namespace qtr {

    SearchData::SearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree,
                           std::shared_ptr<const IdConverter> idConverter, TimeTicker &timeTicker,
                           size_t ansCount,
                           size_t threadCount) :
            ballTree(std::move(ballTree)), idConverter(std::move(idConverter)), timeTicker(timeTicker),
            ansCount(ansCount), threadsCount(threadCount) {}
} // qtr