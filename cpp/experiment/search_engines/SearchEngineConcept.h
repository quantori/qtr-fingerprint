#pragma once

#include <vector>
#include <concepts>

enum class SearchEngineType {
    QtrRDKit,
    QtrIndigo,
    RDKit,
    Indigo,
};

template<typename T>
concept SearchEngine = requires(T t, const std::string &querySmiles, int maxResults, bool &stopFlag) {
    { t.getMatches(querySmiles, maxResults, stopFlag) } -> std::convertible_to<std::vector<uint64_t>>;
};

