#pragma once

#include <concepts>
#include <memory>
#include <vector>
#include <string>

#include "search/utils/SearchQuery.h"
#include "search/utils/SearchResult.h"
#include "frameworks/FrameworkInterface.h"
#include "dataset/SmilesStorage.h"

template<typename SearchEngineT>
concept SearchEngineInterface = requires(SearchEngineT searchEngine) {
    typename SearchEngineT::FrameworkT;
    typename SearchEngineT::ResultT;

    requires FrameworkInterface<typename SearchEngineT::FrameworkT>;

    requires std::constructible_from<SearchEngineT, SmilesStorage &&>;

    {
    searchEngine.search(std::declval<SearchQuery>())
    } -> std::same_as<std::unique_ptr<SearchResult<typename SearchEngineT::ResultT>>>;
};
