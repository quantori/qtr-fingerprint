#pragma once

#include "SearchEngineInterface.h"

#include <memory>

namespace qtr {

struct SearchEngineFactory {
    
    enum SearchEngineType {
        BINGO,
        EXHAUSTIVE,
        DECISION_TREE_TRIVIAL,
        DECISION_TREE_OPTIMAL
    };

    static SearchEnginePtr create(
        SearchEngineType searchEngineType,
        indigo_cpp::IndigoSessionPtr indigoSessionPtr,
        size_t decisionTreeMaxLeafSize = 100);
};

} // namespace qtr