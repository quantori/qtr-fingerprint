#pragma once

#include <filesystem>

#include "base_cpp/array.h"
#include "indigo.h"
#include "frameworks/IndigoFramework.h"

#include "search/engines/SearchEngineInterface.h"


class IndigoSearchEngine {
public:
    using FrameworkT = IndigoFramework;
    using ResultT = SearchResult<FrameworkT>;

    explicit IndigoSearchEngine(SmilesStorage &&dataset);

    [[nodiscard]] std::unique_ptr<ResultT> search(const SearchQuery &query) const;

    ~IndigoSearchEngine();

private:

    IndigoSearchEngine();

    std::filesystem::path _dbFilePath;
    indigo_cpp::BingoMolecule _db;
    std::vector<FrameworkT::MoleculeT> _aliveMolecules; // Trick to make
};

static_assert(SearchEngineInterface<IndigoSearchEngine>, "IndigoSearchEngine must satisfy SearchEngineInterface");
