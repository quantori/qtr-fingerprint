#pragma once

#include "SearchEngine.h"

enum class SearchEngineType {
    QtrRDKit,
    QtrIndigo,
    RDKit,
    Indigo,
};

class SearchEngineFactory {
public:
    static SearchEngine create(SearchEngineType type);
};
