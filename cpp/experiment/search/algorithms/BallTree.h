#pragma once

#include <string>
#include <variant>

#include <glog/logging.h>
#include "GraphMol/SmilesParse/SmilesParse.h"

#include "search/engines/SearchEngineInterface.h"
#include "search/algorithms/FullBinaryTree.h"
#include "frameworks/FrameworkInterface.h"
#include "search/utils/ExtendedSearchQuery.h"
#include "search/utils/utils.h"
#include "dataset/CachedDataset.h"
#include "search/algorithms/BallTreeSplitter.h"
#include "search/utils/BallTreeSearchResult.h"
#include "Profiling.h"

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
    std::optional<NodeData> data;

    BallTreeNode(): centroid(FrameworkT::getEmptyFingerprint()) {}

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
    using QueryMoleculeT = typename FrameworkT::QueryMoleculeT;
    using StorageMoleculeT = typename FrameworkT::StorageMoleculeT;
    using FingerprintT = typename FrameworkT::FingerprintT;
    using Node = BallTreeNode<FrameworkT>;
    using Tree = FullBinaryTree<Node>;
    using ExtendedQueryT = ExtendedSearchQuery<FrameworkT>;
    using ResultT = BallTreeSearchResult<FrameworkT>;

    std::unique_ptr<ResultT > search(const ExtendedQueryT &query) const {
        auto result = std::make_unique<ResultT>();
        size_t nodeId = Tree::root();
        while (nodeId != Tree::endNodeId()) {
            if (checkShouldStopSearch(query, *result)) {
                break;
            }
            bool skipSubtree = shouldSkipSubtree(nodeId, query.fingerprint(), *result);
            if (!skipSubtree && Tree::isLeaf(nodeId)) {
                searchInLeafNode(nodeId, query, *result);
            }
            nodeId = Tree::traverseToNextNode(nodeId, skipSubtree);
        }
        return result;
    }

    explicit BallTree(CachedDataset<FrameworkT> &&dataset, size_t bucketSize) :
            Tree(calculateDepth(dataset.size(), bucketSize)),
            _dataset(std::move(dataset)) {
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

    const CachedDataset<FrameworkT> &dataset() const {
        return _dataset;
    }

private:
    CachedDataset<FrameworkT> _dataset;

    void searchInLeafNode(size_t nodeId, const ExtendedQueryT &query, ResultT &result) const {
        ProfileScope("BallTree searchInLeafNode");
        result.leavesVisited++;
        assert(Tree::isLeaf(nodeId));
        auto &node = this->node(nodeId);
        auto &bucket = node.getLeafData().moleculeIndices;

        for (size_t molIdx: bucket) {
            auto mol = _dataset.molecule(molIdx);
            if (FrameworkT::isSubstructure(query.molecule(), *mol)) {
                result.addResult(*mol);
            }
        }
    }

    bool shouldSkipSubtree(size_t nodeId, const FingerprintT &queryFingerprint, ResultT &result) const {
        if (nodeId == Tree::root() || Tree::isRightChild(nodeId)) {
            return false;
        }
        auto par = Tree::parent(nodeId);
        auto &parNode = this->node(par);
        if (Tree::isLeaf(nodeId)) {
            // check whole fingerprint at leaves
            bool skip = !FrameworkT::isSubFingerprint(queryFingerprint, getCentroid(par));
            result.leavesSKipped += skip;
            return skip;
        }
        // check only a split bit at internal nodes
        size_t splitBit = parNode.getInternalData().splitBit;
        bool skip = FrameworkT::getFingerprintBit(queryFingerprint, splitBit) >
                    FrameworkT::getFingerprintBit(parNode.centroid, splitBit);
        result.internalNodesSkipped += skip;
        return skip;
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
        Node &node = this->node(nodeId);
        assert(!node.isInitialized());
        auto molIndices = splitter.extractIndices(nodeId);
        calculateLeafCentroid(node.centroid, molIndices);
        node.data.emplace(typename Node::LeafData(molIndices));
    }

    void buildInternalNode(size_t nodeId, BallTreeSplitter<FrameworkT> &splitter) {
        assert(!Tree::isLeaf(nodeId));
        Node &node = this->node(nodeId);
        assert(!node.isInitialized());
        calculateInternalNodeCentroid(nodeId);
        size_t splitBit = splitter.getSplitBit(nodeId);
        bool centroidAtSplitBit = FrameworkT::getFingerprintBit(node.centroid, splitBit);
        node.data.emplace(typename Node::InternalData(splitBit, centroidAtSplitBit));
    }

    void calculateLeafCentroid(FingerprintT& centroid, const std::vector<size_t> &moleculeIndices) {
        assert(!moleculeIndices.empty());
        for (size_t j = 0; j < FrameworkT::getFingerprintSize(); j++) {
            for (auto &idx: moleculeIndices) {
                if (FrameworkT::getFingerprintBit(centroid, j)) {
                    break;
                }
                auto &fingerprint = _dataset.fingerprint(idx);
                if (FrameworkT::getFingerprintBit(fingerprint, j)) {
                    FrameworkT::setFingerprintBit(centroid, j, true);
                }
            }
        }
    }

    void calculateInternalNodeCentroid(size_t nodeId) {
        assert(!Tree::isLeaf(nodeId));
        auto& node = this->node(nodeId);
        FingerprintT& centroid = node.centroid;

        auto& fpLeft = getCentroid(Tree::leftChild(nodeId));
        auto &fpRight = getCentroid(Tree::rightChild(nodeId));
        for (size_t idx = 0; idx < FrameworkT::getFingerprintSize(); idx++) {
            if (FrameworkT::getFingerprintBit(fpLeft, idx) || FrameworkT::getFingerprintBit(fpRight, idx)) {
                FrameworkT::setFingerprintBit(centroid, idx, true);
            }
        }
    }

    const FrameworkT::FingerprintT &getCentroid(size_t nodeId) const {
        const Node &node = this->node(nodeId);
        return node.centroid;
    }
};
