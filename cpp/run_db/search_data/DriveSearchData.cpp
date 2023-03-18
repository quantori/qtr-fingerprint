#include "DriveSearchData.h"

#include <utility>

namespace qtr {
    DbType DriveSearchData::getDbType() const {
        return DbType::OnDrive;
    }

    DriveSearchData::DriveSearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree,
                                     std::shared_ptr<const IdConverter> idConverter, TimeTicker &timeTicker,
                                     size_t ansCount, size_t threadCount) :
            SearchData(std::move(ballTree), std::move(idConverter), timeTicker, ansCount, threadCount) {}
} // qtr