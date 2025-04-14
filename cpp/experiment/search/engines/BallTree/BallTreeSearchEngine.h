#pragma once

#include <memory>

#include "search/engines/SearchEngineInterface.h"
#include "frameworks/FrameworkInterface.h"
#include "frameworks/RDKitFramework.h"
#include "search/algorithms/BallTree.h"
#include "search/utils/BallTreeSearchResult.h"

template<typename FrameworkType> requires FrameworkInterface<FrameworkType>
class BallTreeSearchEngine {
public:
    using FrameworkT = FrameworkType;
    using CachedDatasetT = CachedDataset<FrameworkT>;
    using ExtendedSearchQueryT = ExtendedSearchQuery<FrameworkT>;
    using ResultT = BallTreeSearchResult<FrameworkT>;

    static inline const int BucketSize = 2;

    BallTreeSearchEngine() = delete;

    explicit BallTreeSearchEngine(SmilesStorage &&dataset) : _ballTree(CachedDatasetT(std::move(dataset)), BucketSize) {
    }

    [[nodiscard]] std::unique_ptr<ResultT> search(SearchQuery query) const {
        ExtendedSearchQueryT extendedQuery(query);
        auto res = _ballTree.search(extendedQuery);
        return res;
    }

    [[nodiscard]] std::unique_ptr<typename FrameworkT::MoleculeT>
    getMolFromResult(size_t resultIdx, const ResultT &searchResult) const {
        const CachedDataset<FrameworkT> &dataset = _ballTree.dataset();
        size_t datasetIdx = searchResult.get(resultIdx);
        return dataset.molecule(datasetIdx);
    }

private:
    BallTree<FrameworkT> _ballTree;
};

static_assert(SearchEngineInterface<BallTreeSearchEngine<RDKitFramework>>,
              "BallTreeSearchEngine must satisfy SearchEngineInterface");
