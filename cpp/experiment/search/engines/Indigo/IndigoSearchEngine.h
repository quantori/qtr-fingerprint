#pragma once

#include <filesystem>

#include "bingo-nosql.h"

#include "frameworks/IndigoFramework.h"
#include "search/engines/SearchEngineInterface.h"


class IndigoSearchEngine {
public:
    using FrameworkT = IndigoFramework;
    using ResultT = int;

    explicit IndigoSearchEngine(SmilesStorage &&dataset);

    [[nodiscard]] std::unique_ptr<SearchResult<ResultT>> search(const SearchQuery &query) const;

    ~IndigoSearchEngine();

    // TODO: rewrite framework using indigo instead of indigo_cpp? Use custom class to make Molecule copy constructible,
    //      run code with all checks possible
private:

    IndigoSearchEngine();

    std::filesystem::path _dbFilePath;
    indigo_cpp::BingoMolecule _db;
};

static_assert(SearchEngineInterface<IndigoSearchEngine>, "IndigoSearchEngine must satisfy SearchEngineInterface");
