#pragma once

#include <memory>

#include "GraphMol/SubstructLibrary/SubstructLibrary.h"

#include "search/engines/SearchEngineInterface.h"
#include "frameworks/RDKitFramework.h"
#include "dataset/DatasetInterface.h"

class RDKitSearchEngine {
public:
    using FrameworkT = RDKitFramework;
    using ResultT = SearchResult<FrameworkT>;

    RDKitSearchEngine() = delete;

    explicit RDKitSearchEngine(SmilesStorage &&dataset);

    [[nodiscard]] std::unique_ptr<ResultT> search(const SearchQuery &query) const;
private:
    std::unique_ptr<RDKit::SubstructLibrary> _substructLibrary;
};

static_assert(SearchEngineInterface<RDKitSearchEngine>, "RDKitSearchEngine must satisfy SearchEngineInterface");
