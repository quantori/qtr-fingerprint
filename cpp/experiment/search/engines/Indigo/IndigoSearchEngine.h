#pragma once

#include <filesystem>

#include "base_cpp/array.h"
#include "indigo.h"
#include "frameworks/IndigoFramework.h"

#include "search/engines/SearchEngineInterface.h"


class IndigoSearchEngine {
public:
    using FrameworkT = IndigoFramework;

    explicit IndigoSearchEngine(SmilesStorage &&dataset);

    [[nodiscard]] std::unique_ptr<SearchResult> search(const SearchQuery &query) const;

    [[nodiscard]] std::unique_ptr<FrameworkT::MoleculeT>
    getMolFromResult(size_t resultIdx, const SearchResult &searchResult) const;

    ~IndigoSearchEngine();

private:

    IndigoSearchEngine();

    std::filesystem::path _dbFilePath;
    indigo_cpp::BingoMolecule _db;
};

static_assert(SearchEngineInterface<IndigoSearchEngine>, "IndigoSearchEngine must satisfy SearchEngineInterface");
