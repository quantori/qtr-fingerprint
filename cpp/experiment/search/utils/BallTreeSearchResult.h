#pragma once

#include "search/utils/SearchResult.h"

template<typename ResultT>
class BallTreeSearchResult : public SearchResult<ResultT> {
public:
    size_t leavesVisited = 0;
    size_t leavesSKipped = 0;
    size_t internalNodesSkipped = 0;

//    using SearchResult::size;
//    using SearchResult::getNthResult;
//    using SearchResult::addResultByIndex;
};
