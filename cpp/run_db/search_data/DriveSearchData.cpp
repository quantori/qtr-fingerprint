#include "DriveSearchData.h"

#include <utility>

namespace qtr {
    SearchData::DerivedClasses DriveSearchData::getClass() const {
        return SearchData::DerivedClasses::DriveSearchData;
    }

    DriveSearchData::DriveSearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree,
                                     std::shared_ptr<const IdConverter> idConverter, TimeTicker &timeTicker,
                                     size_t ansCount, size_t threadCount, double timeLimit) :
            SearchData(std::move(ballTree), std::move(idConverter), timeTicker, ansCount, threadCount, timeLimit) {}
} // qtr