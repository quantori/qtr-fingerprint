#pragma once

#include "SearchEngine.h"

#include <filesystem>


enum class SearchEngineType {
    QtrRDKit,
    QtrIndigo,
    RDKit,
    Indigo,
};

class SearchEngineFactory {
public:
    static std::unique_ptr<SearchEngine> create(SearchEngineType type, const std::filesystem::path &datasetDir);
};
