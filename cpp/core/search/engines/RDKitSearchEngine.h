#pragma once

#include <memory>

#include "GraphMol/SubstructLibrary/SubstructLibrary.h"

#include "search/engines/SearchEngineInterface.h"
#include "frameworks/RDKitFramework.h"
#include "dataset/DatasetInterface.h"
#include "utils/Config.h"

class RDKitSearchEngine {
public:
    using FrameworkT = RDKitFramework;
    using ResultT = size_t;

    RDKitSearchEngine() = delete;

    explicit RDKitSearchEngine(FrameworkT framework, SmilesStorage &&dataset, const Config &config);

    [[nodiscard]] std::unique_ptr<SearchResult<ResultT>> search(const SearchQuery &query) const;

    static StatTable getStat();

    [[nodiscard]] std::string resultToSmiles(const ResultT &result) const;

private:
    FrameworkT _framework;
    std::unique_ptr<RDKit::SubstructLibrary> _substructLibrary;
};

static_assert(SearchEngineInterface<RDKitSearchEngine>, "RDKitSearchEngine must satisfy SearchEngineInterface");
