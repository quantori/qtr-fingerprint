#pragma once

#include "search/utils/SearchQuery.h"
#include "search/utils/SearchResult.h"

template<typename ResultT>
inline bool checkShouldStopSearch(const SearchQuery &query, const SearchResult<ResultT> &result) {
    return query.checkStopFlag() || result.size() >= query.maxResults();
}

