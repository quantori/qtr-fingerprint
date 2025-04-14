#pragma once

#include "search/utils/SearchQuery.h"
#include "search/utils/SearchResult.h"

template<typename FrameworkT>
inline bool checkShouldStopSearch(const SearchQuery& query, const SearchResult<FrameworkT>& result) {
    return query.checkStopFlag() || result.size() >= query.maxResults();
}