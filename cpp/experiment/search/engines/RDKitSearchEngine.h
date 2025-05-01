#pragma once

#include <memory>

#include "GraphMol/SubstructLibrary/SubstructLibrary.h"

#include "search/engines/SearchEngineInterface.h"
#include "frameworks/RDKitFramework.h"
#include "dataset/DatasetInterface.h"
#include "search/utils/SearchEngineConfig.h"

class RDKitSearchEngine {
public:
    using FrameworkT = RDKitFramework;
    using ResultT = FrameworkT::MoleculeT;

    RDKitSearchEngine() = delete;

    explicit RDKitSearchEngine(SmilesStorage &&dataset, const SearchEngineConfig& config);

    [[nodiscard]] std::unique_ptr<SearchResult<ResultT>> search(const SearchQuery &query) const;

    StatTable getStat() const;
private:
    std::unique_ptr<RDKit::SubstructLibrary> _substructLibrary;
};

static_assert(SearchEngineInterface<RDKitSearchEngine>, "RDKitSearchEngine must satisfy SearchEngineInterface");
