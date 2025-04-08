#pragma once

#include "search/utils/SearchResult.h"

class BallTreeSearchResult : public SearchResult {
public:
    size_t leavesVisited = 0;
    size_t leavesSKipped = 0;
    size_t internalNodesSkipped = 0;

//    using SearchResult::size;
//    using SearchResult::getNthResult;
//    using SearchResult::addResultByIndex;
};
