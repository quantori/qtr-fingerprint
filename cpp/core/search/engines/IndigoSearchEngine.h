#pragma once

#include <filesystem>

#include "bingo-nosql.h"

#include "frameworks/IndigoCppFramework.h"
#include "search/engines/SearchEngineInterface.h"
#include "utils/Config.h"


class IndigoSearchEngine {
public:
    using FrameworkT = IndigoCppFramework;
    using ResultT = size_t;

    explicit IndigoSearchEngine(FrameworkT framework, SmilesStorage &&dataset, const Config &config);

    [[nodiscard]] std::unique_ptr<SearchResult<ResultT>> search(const SearchQuery &query) const;

    ~IndigoSearchEngine();

    [[nodiscard]] StatTable getStat() const;

    std::string resultToSmiles(const ResultT &result) const;

private:

    explicit IndigoSearchEngine(FrameworkT framework);

    FrameworkT _framework;
    std::filesystem::path _dbFilePath;
    indigo_cpp::BingoMolecule _db;
};

static_assert(SearchEngineInterface<IndigoSearchEngine>, "IndigoSearchEngine must satisfy SearchEngineInterface");
