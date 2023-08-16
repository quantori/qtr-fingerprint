#pragma once

#include <memory>

#include "SearchData.h"
#include "RunArgs.h"
#include "TimeTicker.h"

namespace qtr {

    class SearchDataLoader {
    public:
        static std::shared_ptr<SearchData> load(const RunArgs &args, TimeTicker &timeTicker);
    };

} // qtr 
