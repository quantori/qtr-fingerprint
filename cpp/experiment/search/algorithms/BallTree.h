#pragma once

#include <string>
#include <variant>

#include <glog/logging.h>
#include "GraphMol/SmilesParse/SmilesParse.h"

#include "search/engines/SearchEngineInterface.h"
#include "search/algorithms/FullBinaryTree.h"
#include "frameworks/FrameworkInterface.h"
#include "search/utils/ExtendedSearchQuery.h"
#include "search/utils/SearchUtils.h"
#include "dataset/CachedDataset.h"
#include "search/algorithms/BallTreeSplitter.h"
#include "Profiling.h"
#include "BallTreeQueryStat.h"
#include "BallTreeFingerprintChecker.h"
#include "BallTreeNodesStat.h"

template<typename FrameworkT> requires FrameworkInterface<FrameworkT>
class BallTreeNode {
public:
    using FingerprintT = FrameworkT::FingerprintT;

    struct LeafData {
        std::vector<size_t> moleculeIndices;
    };

    struct InternalData {
        size_t splitBit;
        bool centroidAtSplitBit;
    };

    using NodeData = std::variant<LeafData, InternalData>;

    FingerprintT centroid;
    size_t subsetSize;
    std::optional<NodeData> data;

    BallTreeNode() : centroid(FrameworkT::getEmptyFingerprint()) {}

    [[nodiscard]] bool isLeaf() const {
        return data && std::holds_alternative<LeafData>(*data);
    }

    [[nodiscard]] bool isInitialized() const {
        return data.has_value();
    }

    const LeafData &getLeafData() const {
        return std::get<LeafData>(*data);
    }

    const InternalData &getInternalData() const {
        return std::get<InternalData>(*data);
    }

    explicit BallTreeNode(FingerprintT &&centroid)
            : centroid(std::move(centroid)) {}
};

template<typename FrameworkT> requires FrameworkInterface<FrameworkT>
class BallTree : public FullBinaryTree<BallTreeNode<FrameworkT>> {
public:
    using MoleculeT = typename FrameworkT::MoleculeT;
    using FingerprintT = typename FrameworkT::FingerprintT;
    using NodeT = BallTreeNode<FrameworkT>;
    using Tree = FullBinaryTree<NodeT>;
    using ExtendedQueryT = ExtendedSearchQuery<FrameworkT>;
    using ResultT = size_t;

    using Tree::rightChild;
    using Tree::leftChild;
    using Tree::node;
    using Tree::root;
    using Tree::endNodeId;
    using Tree::traverseToNextNode;
    using Tree::isLeaf;

    std::unique_ptr<SearchResult<ResultT>> search(const ExtendedQueryT &query) const {
        ProfileScope("BallTree::search");
        auto result = std::make_unique<SearchResult<ResultT>>();
        BallTreeFingerprintChecker<FrameworkT> fpChecker(query.fingerprint());
        BallTreeQueryStat stat(Tree::depth());
        size_t nodeId = Tree::root();
        while (nodeId != Tree::endNodeId()) {
            if (checkShouldStopSearch(query, *result)) {
                break;
            }
            bool skipSubtree = shouldSkipSubtree(nodeId, fpChecker, *result, stat);
            if (!skipSubtree) {
                _nodesStat[nodeId].visits++;
                stat.nodesVisitedPerDepth.at(Tree::nodeDepth(nodeId))++;
                stat.subsetSizePerDepth.at(Tree::nodeDepth(nodeId)) += this->node(nodeId).subsetSize;
                if (Tree::isLeaf(nodeId)) {
                    searchInLeafNode(nodeId, query, *result, stat);
                }
            } else {
                _nodesStat[nodeId].skips++;
            }
            nodeId = Tree::traverseToNextNode(nodeId, skipSubtree);
        }
        result->setStatRow(stat.toStatRow());
        return result;
    }

    explicit BallTree(CachedDataset<FrameworkT> &&dataset, size_t bucketSize) :
            Tree(calculateDepth(dataset.size(), bucketSize)),
            _dataset(std::move(dataset)), _nodesStat(Tree::nodeCount()) {
        BallTreeSplitter<FrameworkT> splitter(_dataset, Tree::root(), Tree::nodeCount());
        for (size_t nodeId = Tree::root(); nodeId < Tree::nodeCount(); nodeId++) {
            if (Tree::isLeaf(nodeId)) {
                buildLeaf(nodeId, splitter);
            } else {
                splitter.split(nodeId, Tree::leftChild(nodeId), Tree::rightChild(nodeId));
            }
        }
        for (int nodeId = (int) Tree::lastInternalNodeId(); nodeId >= (int) Tree::root(); nodeId--) {
            buildInternalNode(nodeId, splitter);
        }
    }

    StatTable getNodesStat() const {
        return _nodesStat.toStatTable();
    }

    const CachedDataset<FrameworkT> &dataset() const {
        return _dataset;
    }

private:
    CachedDataset<FrameworkT> _dataset;
    mutable BallTreeNodesStat _nodesStat;

    void searchInLeafNode(size_t nodeId, const ExtendedQueryT &query, SearchResult<ResultT> &result,
                          BallTreeQueryStat &queryStat) const {
        ProfileScope("BallTree::searchInLeafNode");
        queryStat.leafSearches++;
        assert(Tree::isLeaf(nodeId));
        auto &node = this->node(nodeId);
        auto &bucket = node.getLeafData().moleculeIndices;

        for (size_t molIdx: bucket) {
            auto &fp = _dataset.fingerprint(molIdx);
            if (!FrameworkT::isSubFingerprint(query.fingerprint(), fp)) {
                continue;
            }
            auto mol = _dataset.molecule(molIdx);
            if (!FrameworkT::isSubstructure(query.molecule(), *mol)) {
                continue;
            }
            result.addResult(molIdx);
        }
    }

    bool shouldSkipSubtree(size_t nodeId, const BallTreeFingerprintChecker<FrameworkT> &fpChecker,
                           SearchResult<ResultT> &result, BallTreeQueryStat &stat) const {
        if (nodeId == Tree::root() || Tree::isRightChild(nodeId)) {
            return false;
        }
        return !fpChecker.check(getCentroid(nodeId));
    }

    static size_t calculateDepth(size_t datasetSize, size_t bucketSize) {
        auto depth = (size_t) std::ceil(std::log2(datasetSize / bucketSize));
        if (depth < 2) {
            throw std::runtime_error("Cannot build BallTree with depth < 2");
        }
        LOG(INFO) << "BallTree selected depth = " << depth << " with >= " << datasetSize / (1ull << depth)
                  << " items in bucket";
        return depth;
    }

    void buildLeaf(size_t nodeId, BallTreeSplitter<FrameworkT> &splitter) {
        assert(Tree::isLeaf(nodeId));
        auto &node = this->node(nodeId);
        assert(!node.isInitialized());
        auto molIndices = splitter.extractIndices(nodeId);
        node.subsetSize = molIndices.size();
        calculateLeafCentroid(node.centroid, molIndices);
        node.data.emplace(typename NodeT::LeafData(molIndices));
    }

    void buildInternalNode(size_t nodeId, BallTreeSplitter<FrameworkT> &splitter) {
        assert(!Tree::isLeaf(nodeId));
        auto &node = this->node(nodeId);
        assert(!node.isInitialized());
        calculateInternalNodeCentroid(nodeId);
        size_t splitBit = splitter.getSplitBit(nodeId);
        bool centroidAtSplitBit = FrameworkT::getFingerprintBit(node.centroid, splitBit);
        node.subsetSize =
                this->node(Tree::leftChild(nodeId)).subsetSize + this->node(Tree::rightChild(nodeId)).subsetSize;
        node.data.emplace(typename NodeT::InternalData(splitBit, centroidAtSplitBit));
    }

    void calculateLeafCentroid(FingerprintT &centroid, const std::vector<size_t> &moleculeIndices) {
        assert(!moleculeIndices.empty());
        for (size_t j = 0; j < FrameworkT::getFingerprintSize(); j++) {
            for (auto &idx: moleculeIndices) {
                auto &fingerprint = _dataset.fingerprint(idx);
                if (FrameworkT::getFingerprintBit(fingerprint, j)) {
                    FrameworkT::setFingerprintBit(centroid, j, true);
                    break;
                }
            }
        }
    }

    void calculateInternalNodeCentroid(size_t nodeId) {
        assert(!Tree::isLeaf(nodeId));
        auto &node = this->node(nodeId);
        FingerprintT &centroid = node.centroid;

        auto &fpLeft = getCentroid(Tree::leftChild(nodeId));
        auto &fpRight = getCentroid(Tree::rightChild(nodeId));
        for (size_t idx = 0; idx < FrameworkT::getFingerprintSize(); idx++) {
            if (FrameworkT::getFingerprintBit(fpLeft, idx) || FrameworkT::getFingerprintBit(fpRight, idx)) {
                FrameworkT::setFingerprintBit(centroid, idx, true);
            }
        }
    }

    const FrameworkT::FingerprintT &getCentroid(size_t nodeId) const {
        const auto &node = this->node(nodeId);
        return node.centroid;
    }
};
