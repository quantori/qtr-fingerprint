#pragma once

#include <string>

#include "answer_filtering/PropertiesFilter.h"
#include "BallTreeDriveSearchEngine.h"
#include "modes/web/IdConverter.h"
#include "DbConfig.h"

namespace qtr {

    class SearchData {
    public:

        [[nodiscard]] virtual DbType getDbType() const = 0;

        virtual ~SearchData() = default;

        SearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree, std::shared_ptr<const IdConverter> idConverter,
                   TimeTicker &timeTicker, size_t ansCount,
                   size_t threadCount);


        std::shared_ptr<const BallTreeSearchEngine> ballTree;
        std::shared_ptr<const IdConverter> idConverter;
        TimeTicker &timeTicker;
        size_t ansCount;
        size_t threadsCount;
    };

} // qtr

