#pragma once

#include <string>

#include "answer_filtering/PropertiesFilter.h"
#include "BallTreeDriveSearchEngine.h"
#include "modes/web/IdConverter.h"

namespace qtr {

    class SearchData {
    public:
        virtual ~SearchData() = default;

        SearchData(size_t ansCount, size_t threadCount, double timeLimit);

        virtual std::unique_ptr<QueryData<CIDType>>
        search(const std::string &querySmiles, const PropertiesFilter::Bounds &queryBounds) = 0;

        size_t ansCount;
        size_t threadsCount;
        double timeLimit;
    };

} // qtr