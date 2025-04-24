#pragma once

#include <memory>

#include "search/engines/SearchEngineInterface.h"
#include "frameworks/FrameworkInterface.h"
#include "frameworks/RDKitFramework.h"
#include "search/algorithms/BallTree.h"

template<typename FrameworkType> requires FrameworkInterface<FrameworkType>
class BallTreeSearchEngine {
public:
    using FrameworkT = FrameworkType;
    using CachedDatasetT = CachedDataset<FrameworkT>;
    using ExtendedSearchQueryT = ExtendedSearchQuery<FrameworkT>;
    using ResultT = BallTree<FrameworkT>::ResultT;

    static inline const int BucketSize = 3;

    BallTreeSearchEngine() = delete;

    explicit BallTreeSearchEngine(SmilesStorage &&dataset) : _ballTree(CachedDatasetT(std::move(dataset)), BucketSize) {
    }

    [[nodiscard]] std::unique_ptr<SearchResult<ResultT>> search(SearchQuery query) const {
        ExtendedSearchQueryT extendedQuery(query);
        auto res = _ballTree.search(extendedQuery);
        return res;
    }

    [[nodiscard]] StatRow getStat() const;

private:
    BallTree<FrameworkT> _ballTree;
};

template<typename FrameworkType>
requires FrameworkInterface<FrameworkType>StatRow BallTreeSearchEngine<FrameworkType>::getStat() const {
    // TODO
    return StatRow();
}

static_assert(SearchEngineInterface<BallTreeSearchEngine<RDKitFramework>>,
              "BallTreeSearchEngine must satisfy SearchEngineInterface");
