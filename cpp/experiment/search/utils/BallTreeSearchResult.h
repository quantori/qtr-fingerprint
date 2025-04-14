#pragma once

#include "search/utils/SearchResult.h"

template<typename FrameworkT>
class BallTreeSearchResult : public SearchResult<FrameworkT> {
public:
    size_t leavesVisited = 0;
    size_t leavesSKipped = 0;
    size_t internalNodesSkipped = 0;

//    using SearchResult::size;
//    using SearchResult::getNthResult;
//    using SearchResult::addResultByIndex;
};
