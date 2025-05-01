#pragma once

#include <filesystem>

#include "bingo-nosql.h"

#include "frameworks/IndigoCppFramework.h"
#include "search/engines/SearchEngineInterface.h"
#include "search/utils/SearchEngineConfig.h"


class IndigoSearchEngine {
public:
    using FrameworkT = IndigoCppFramework;
    using ResultT = int;

    explicit IndigoSearchEngine(SmilesStorage &&dataset, const SearchEngineConfig& config);

    [[nodiscard]] std::unique_ptr<SearchResult<ResultT>> search(const SearchQuery &query) const;

    ~IndigoSearchEngine();

    StatTable getStat() const;
private:

    IndigoSearchEngine();

    std::filesystem::path _dbFilePath;
    indigo_cpp::BingoMolecule _db;
};

static_assert(SearchEngineInterface<IndigoSearchEngine>, "IndigoSearchEngine must satisfy SearchEngineInterface");
