#pragma once

#include <vector>
#include <cassert>

template<typename NodeType>
class FullBinaryTree {
protected:
    size_t _depth;
    std::vector<NodeType> _nodes;

    // Tree structure methods
    [[nodiscard]] static size_t leftChild(size_t nodeId) {
        return nodeId * 2 + 1;
    }

    [[nodiscard]] static size_t rightChild(size_t nodeId) {
        return nodeId * 2 + 2;
    }

    [[nodiscard]] static size_t parent(size_t nodeId) {
        return (nodeId - 1) >> 1;
    }

    [[nodiscard]] static bool isLeftChild(size_t nodeId) {
        return nodeId != root() && (nodeId & 1) == 1;
    }

    [[nodiscard]] static bool isRightChild(size_t nodeId) {
        return nodeId != root() && (nodeId & 1) == 0;
    }

    [[nodiscard]] static size_t root() {
        return 0;
    }

    [[nodiscard]] static size_t endNodeId() {
        return -1;
    }

    [[nodiscard]] size_t lastInternalNodeId() const {
        return (1u << _depth) - 2;
    }

    [[nodiscard]] static size_t nodeDepth(size_t nodeId) {
        int clz = __builtin_clzll(static_cast<unsigned long long>(nodeId) + 1ull);
        assert(clz >= 0 && clz <= 63);
        return 63 - clz;
    }

    [[nodiscard]] bool isLeaf(size_t nodeId) const {
        return leftChild(nodeId) >= _nodes.size();
    }

    [[nodiscard]] size_t traverseUp(size_t nodeId) const {
        while (true) {
            if (nodeId == root()) {
                return endNodeId();
            }
            size_t par = parent(nodeId);
            if (nodeId == leftChild(par)) {
                return rightChild(par);
            } else {
                nodeId = par;
            }
        }
    }

    [[nodiscard]] size_t traverseDown(size_t nodeId) const {
        if (isLeaf(nodeId)) {
            return traverseUp(nodeId);
        } else {
            return leftChild(nodeId);
        }
    }

    [[nodiscard]] size_t traverseToNextNode(size_t nodeId, bool skipSubtree = false) const {
        if (skipSubtree) {
            return traverseUp(nodeId);
        } else {
            return traverseDown(nodeId);
        }
    }

    [[nodiscard]] size_t nodeIdToLeafId(size_t leafId) const {
        assert((1ull << _depth) - 1 <= leafId);
        leafId -= (1ull << _depth) - 1;
        assert(leafId < (1ull << _depth));
        return leafId;
    }

public:
    explicit FullBinaryTree(size_t depth)
            : _depth(depth), _nodes((2ull << depth) - 1) {
    }

    [[nodiscard]] size_t depth() const { return _depth; }

    [[nodiscard]] size_t nodeCount() const { return _nodes.size(); }

    [[nodiscard]] size_t leafCount() const { return 1ull << _depth; }

    [[nodiscard]] const NodeType &node(size_t nodeId) const {
        assert(nodeId < _nodes.size());
        return _nodes[nodeId];
    }

    [[nodiscard]] NodeType &node(size_t nodeId) {
        assert(nodeId < _nodes.size());
        return _nodes[nodeId];
    }
};
