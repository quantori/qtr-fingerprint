#pragma once

#include <memory>

#include "SearchData.h"
#include "RunArgs.h"

namespace qtr {

    class SearchDataLoader {
    public:
        static std::shared_ptr<SearchData> load(const RunArgs &args);
    };

} // qtr 
a