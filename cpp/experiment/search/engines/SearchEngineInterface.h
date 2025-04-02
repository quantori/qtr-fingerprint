#pragma once

#include <concepts>
#include <memory>
#include <vector>
#include <string>

#include "search/utils/SearchQuery.h"
#include "search/utils/SearchResult.h"
#include "frameworks/FrameworkInterface.h"
#include "dataset/SmilesStorage.h"

template<typename T>
concept DerivedSearchResultPtr =
requires {
    typename T::element_type;
} &&
std::derived_from<typename T::element_type, SearchResult> &&
std::same_as<T, std::unique_ptr<typename T::element_type>>;

template<typename SearchEngineT>
concept SearchEngineInterface = requires(SearchEngineT searchEngine) {
    typename SearchEngineT::FrameworkT;

    requires FrameworkInterface<typename SearchEngineT::FrameworkT>;
    requires std::constructible_from<SearchEngineT, SmilesStorage &&>;

    {
    searchEngine.search(std::declval<SearchQuery>())
    } -> DerivedSearchResultPtr;

    {
    searchEngine.getMolFromResult(std::declval<size_t>(), std::declval<const SearchResult &>())
    } -> std::same_as<std::unique_ptr<typename SearchEngineT::FrameworkT::MoleculeT>>;
};
