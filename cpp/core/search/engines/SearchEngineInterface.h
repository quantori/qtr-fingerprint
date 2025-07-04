#pragma once

#include <concepts>
#include <memory>
#include <vector>
#include <string>

#include "search/utils/SearchQuery.h"
#include "search/utils/SearchResult.h"
#include "utils/Config.h"
#include "frameworks/FrameworkInterface.h"
#include "dataset/SmilesStorage.h"
#include "stats/StatRow.h"

template<typename SearchEngineT>
concept SearchEngineInterface = requires(SearchEngineT searchEngine) {
    typename SearchEngineT::FrameworkT;
    typename SearchEngineT::ResultT;

    requires FrameworkInterface<typename SearchEngineT::FrameworkT>;

    requires std::constructible_from<SearchEngineT, typename SearchEngineT::FrameworkT, SmilesStorage &&, const Config &>;

    {
    searchEngine.search(std::declval<SearchQuery>())
    } -> std::same_as<std::unique_ptr<SearchResult<typename SearchEngineT::ResultT>>>;

    {
    searchEngine.getStat()
    } -> std::same_as<StatTable>;

    {
    searchEngine.resultToSmiles(std::declval<const typename SearchEngineT::ResultT &>())
    } -> std::same_as<std::string>;
};
