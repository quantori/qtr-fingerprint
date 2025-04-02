#pragma once

#include "search/utils/SearchQuery.h"
#include "search/utils/SearchResult.h"

inline bool checkShouldStopSearch(const SearchQuery& query, const SearchResult& result) {
    return query.checkStopFlag() || result.size() >= query.maxResults();
}