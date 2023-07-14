#pragma once

#include <string>

#include "answer_filtering/PropertiesFilter.h"
#include "BallTreeDriveSearchEngine.h"
#include "modes/web/IdConverter.h"
#include "TimeTicker.h"

namespace qtr {

    class SearchData {
    public:

        enum class DerivedClasses {
            RamSearchData,
            DriveSearchData
        };

        [[nodiscard]] virtual DerivedClasses getClass() const = 0;

        virtual ~SearchData() = default;

        SearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree, std::shared_ptr<const IdConverter> idConverter,
                   TimeTicker &timeTicker, size_t ansCount, size_t threadCount, double timeLimit);


        std::shared_ptr<const BallTreeSearchEngine> ballTree;
        std::shared_ptr<const IdConverter> idConverter;
        TimeTicker &timeTicker;
        size_t ansCount;
        size_t threadsCount;
        double timeLimit;
    };

} // qtr

