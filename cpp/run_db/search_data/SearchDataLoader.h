#pragma once

#include <memory>

#include "SearchData.h"
#include "Args.h"
#include "TimeTicker.h"

namespace qtr {

    class SearchDataLoader {
    public:
        static std::shared_ptr<SearchData> load(const Args &args, TimeTicker &timeTicker);
    };

} // qtr 

