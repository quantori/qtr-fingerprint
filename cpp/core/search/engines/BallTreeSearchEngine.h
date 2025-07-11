#pragma once

#include <memory>

#include "search/engines/SearchEngineInterface.h"
#include "frameworks/FrameworkInterface.h"
#include "frameworks/RDKitFramework.h"
#include "search/algorithms/BallTree.h"
#include "utils/Config.h"

template<typename FrameworkType> requires FrameworkInterface<FrameworkType>
class BallTreeSearchEngine {
public:
    using FrameworkT = FrameworkType;
    using CachedDatasetT = CachedDataset<FrameworkT>;
    using ExtendedSearchQueryT = ExtendedSearchQuery<FrameworkT>;
    using ResultT = BallTree<FrameworkT>::ResultT;

    static inline const int BucketSize = 4;

    BallTreeSearchEngine() = delete;

    explicit BallTreeSearchEngine(FrameworkT framework, SmilesStorage &&dataset, const Config &config)
            : _framework(framework),
              _ballTree(CachedDatasetT(std::move(dataset)), BucketSize, config.getInt("depth", -1)) {
    }

    [[nodiscard]] std::unique_ptr<SearchResult<ResultT>> search(SearchQuery query) const {
        ExtendedSearchQueryT extendedQuery(query);
        auto res = _ballTree.search(extendedQuery);
        return res;
    }

    [[nodiscard]] StatTable getStat() const {
        StatTable statTable = _ballTree.getNodesStat();
        finalizedNodesStat(statTable);
        return statTable;
    }

    [[nodiscard]] std::string resultToSmiles(const ResultT &result) const {
        const auto &dataset = _ballTree.dataset();
        auto mol = dataset.molecule(result);
        auto smiles = _framework.moleculeToSmiles(*mol);
        return smiles;
    }

    void finalizedNodesStat(StatTable &statTable) const {
        const auto &dataset = _ballTree.dataset();
        for (size_t nodeId = _ballTree.root();
             nodeId != _ballTree.endNodeId(); nodeId = _ballTree.traverseToNextNode(nodeId)) {

            const typename BallTree<FrameworkT>::NodeT &node = _ballTree.node(nodeId);
            std::vector<const typename FrameworkT::FingerprintT *> fingerprints;
            if (_ballTree.isLeaf(nodeId)) {
                for (auto &molIdx: node.getLeafData().moleculeIndices) {
                    fingerprints.emplace_back(&dataset.fingerprint(molIdx));
                }
            } else {
                const typename BallTree<FrameworkT>::NodeT &lChild = _ballTree.node(_ballTree.leftChild(nodeId));
                const typename BallTree<FrameworkT>::NodeT &rChild = _ballTree.node(_ballTree.rightChild(nodeId));
                fingerprints.emplace_back(&lChild.centroid);
                fingerprints.emplace_back(&rChild.centroid);
            }
            size_t diffBits = countDiffBits<FrameworkT>(fingerprints);
            statTable[nodeId].addEntry("childrenFingerprintDiffBits", diffBits);
        }
    }

    template<typename FrameworkT>
    requires FrameworkInterface<FrameworkT>
    static size_t countDiffBits(std::vector<const typename FrameworkT::FingerprintT *> &fingerprints) {
        size_t result = 0;
        auto &framework = FrameworkT::getInstance();
        for (size_t i = 0; i < framework.getFingerprintSize(); i++) {
            bool hasZero = false;
            bool hasOne = false;
            for (size_t j = 0; j < fingerprints.size() && (!hasZero || !hasOne); j++) {
                if (FrameworkT::getFingerprintBit(*fingerprints[j], i)) {
                    hasOne = true;
                } else {
                    hasZero = true;
                }
            }
            result += hasOne && hasZero;
        }
        return result;
    }

private:
    FrameworkT _framework;
    BallTree<FrameworkT> _ballTree;
};


static_assert(SearchEngineInterface<BallTreeSearchEngine<RDKitFramework>>,
              "BallTreeSearchEngine must satisfy SearchEngineInterface");
