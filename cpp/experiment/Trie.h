#pragma once

#include "IndigoUtils.h"

template<typename MoleculeType, typename QueryMoleculeType, typename FingerprintType>
class Trie {
private:
    struct Node {
        bool isLeaf;
        std::vector<size_t> oneBitsPositions;

        union NodeData {
            struct InnerNode {
                std::unique_ptr<Node> left;
                std::unique_ptr<Node> right;
                size_t bitIndex;
            } inner;

            struct LeafNode {
                std::vector<std::pair<std::unique_ptr<MoleculeType>, std::unique_ptr<FingerprintType>>> molecules;
            } leaf;

            NodeData() {}

            ~NodeData() {}
        } data;

        explicit Node(bool isLeaf) : isLeaf(isLeaf) {
            if (isLeaf) {
                new(&data.leaf) NodeData::LeafNode();
            } else {
                new(&data.inner) NodeData::InnerNode();
            }
        }

        ~Node() {
            if (isLeaf) {
                data.leaf.molecules.~vector();
            } else {
                data.inner.left.~unique_ptr();
                data.inner.right.~unique_ptr();
            }
        }
    };

    std::unique_ptr<Node> _root;
    size_t _maxDepth;

    void buildTrie(Node &node,
                   std::vector<std::pair<std::unique_ptr<MoleculeType>, std::unique_ptr<FingerprintType>>> &data,
                   size_t leftSplitPos = 0,
                   size_t depth = 0) {
        assert(!data.empty());

        const size_t fpSize = data.front().second->size();

        auto [rightSplitPos, zeroCount, oneBitsPositions] = findSplitBit(data, leftSplitPos, fpSize);

        if (rightSplitPos == fpSize || depth >= _maxDepth) {
            createLeafNode(node, std::move(data), std::move(oneBitsPositions));
            return;
        }

        splitData(node, data, rightSplitPos, depth, zeroCount);
    }

    std::tuple<size_t, size_t, std::vector<size_t>>
    findSplitBit(const std::vector<std::pair<std::unique_ptr<MoleculeType>, std::unique_ptr<FingerprintType>>> &data,
                 size_t splitBit, size_t fpSize) {
        std::vector<size_t> oneBitsPositions;
        size_t zeroCount = 0;
        for (; splitBit < fpSize; ++splitBit) {
            zeroCount = countZeroBits(data, splitBit);

            if (zeroCount > 0 && zeroCount < data.size()) {
                break;
            }

            if (zeroCount == data.size()) {
                oneBitsPositions.push_back(splitBit);
            }
        }
        return {splitBit, zeroCount, oneBitsPositions};
    }

    size_t
    countZeroBits(const std::vector<std::pair<std::unique_ptr<MoleculeType>, std::unique_ptr<FingerprintType>>> &data,
                  size_t splitBit) {
        size_t zeroCount = 0;
        for (const auto &[mol, fp]: data) {
            zeroCount += fp->getBit(splitBit) == 0;
        }
        return zeroCount;
    }

    void createLeafNode(Node &node,
                        std::vector<std::pair<std::unique_ptr<MoleculeType>, std::unique_ptr<FingerprintType>>> &&data,
                        std::vector<size_t> &&oneBitsPositions) {
        node.isLeaf = true;
        new(&node.data.leaf) Node::NodeData::LeafNode();
        node.data.leaf.molecules = std::move(data);
        node.oneBitsPositions = std::move(oneBitsPositions);
    }

    void splitData(Node &node,
                   std::vector<std::pair<std::unique_ptr<MoleculeType>, std::unique_ptr<FingerprintType>>> &data,
                   size_t splitBit, size_t depth, size_t zeroCount) {
        auto [zeroData, oneData] = splitDataByBit(data, splitBit, zeroCount);

        node.isLeaf = false;
        node.data.inner.bitIndex = splitBit;
        node.data.inner.left = std::make_unique<Node>(true);
        node.data.inner.right = std::make_unique<Node>(true);

        // Build oneBitsPositions for the inner node
        for (size_t i = 0; i < splitBit; ++i) {
            if (data.front().second->getBit(i) == 1) {
                node.oneBitsPositions.push_back(i);
            }
        }

        buildTrie(*node.data.inner.left, zeroData, splitBit + 1, depth + 1);
        buildTrie(*node.data.inner.right, oneData, splitBit + 1, depth + 1);
    }

    std::pair<std::vector<std::pair<std::unique_ptr<MoleculeType>, std::unique_ptr<FingerprintType>>>,
            std::vector<std::pair<std::unique_ptr<MoleculeType>, std::unique_ptr<FingerprintType>>>>
    splitDataByBit(std::vector<std::pair<std::unique_ptr<MoleculeType>, std::unique_ptr<FingerprintType>>> &data,
                   size_t splitBit, size_t zeroCount) {
        std::vector<std::pair<std::unique_ptr<MoleculeType>, std::unique_ptr<FingerprintType>>> zeroData;
        std::vector<std::pair<std::unique_ptr<MoleculeType>, std::unique_ptr<FingerprintType>>> oneData;
        zeroData.reserve(zeroCount);
        oneData.reserve(data.size() - zeroCount);
        for (auto &[mol, fp]: data) {
            if (fp->getBit(splitBit) == 0) {
                zeroData.emplace_back(std::move(mol), std::move(fp));
            } else {
                oneData.emplace_back(std::move(mol), std::move(fp));
            }
        }
        assert(zeroData.capacity() == zeroData.size());
        assert(oneData.capacity() == oneData.size());

        data.clear();
        data.shrink_to_fit();

        return {std::move(zeroData), std::move(oneData)};
    }

    void search(QueryMoleculeType &mol, const FingerprintType &fingerprint, int maxResults, bool &stopFlag,
                Node &node, std::vector<uint64_t> &results) const {
        std::vector<std::pair<Node *, size_t>> stack;
        stack.emplace_back(&node, 0);

        while (!stack.empty() && (int64_t)results.size() < maxResults && !stopFlag) {
            auto [current, lastSplitBit] = stack.back();
            stack.pop_back();

            if (canSkipNode(*current, fingerprint, lastSplitBit)) {
                continue;
            }

            if (current->isLeaf) {
                processLeafNode(*current, mol, maxResults, stopFlag, results);
            } else {
                processInnerNode(*current, fingerprint, lastSplitBit, stack);
            }
        }
    }

    bool canSkipNode(const Node &node, const FingerprintType &fingerprint, size_t startBitPos) const {
        auto &oneBitsPositions = node.oneBitsPositions;
        size_t lastBitPos = node.isLeaf ? fingerprint.size() : node.data.inner.bitIndex;
        auto oneBitIterator = oneBitsPositions.cbegin();
        for (size_t bitPos = startBitPos; bitPos < lastBitPos && oneBitIterator != oneBitsPositions.end(); bitPos++) {
            bool nodeVal = *oneBitIterator == bitPos;
            bool fpVal = fingerprint.getBit(bitPos) == 1;
            if (fpVal < nodeVal) {
                return true;
            }
            if (nodeVal) {
                oneBitIterator++;
            }
        }
        return false;
    }

    void processLeafNode(const Node &node, QueryMoleculeType &mol, int maxResults, bool &stopFlag,
                         std::vector<uint64_t> &results) const {
        for (const auto &[storedMol, storedFp]: node.data.leaf.molecules) {
            if ((int64_t) results.size() >= maxResults || stopFlag) {
                return;
            }
            if (isSubstructure(mol, *storedMol)) {
                // TODO: use indices or other information about molecules
                results.push_back(-2); // -2 is a random selected number
            }
        }
    }

    void processInnerNode(const Node &node, const FingerprintType &fingerprint, size_t lastBitPos,
                          std::vector<std::pair<Node *, size_t>> &stack) const {
        stack.emplace_back(node.data.inner.right.get(), lastBitPos);
        if (fingerprint.getBit(node.data.inner.bitIndex) == 0) {
            stack.emplace_back(node.data.inner.left.get(), lastBitPos);
        }
    }

public:
    std::vector<uint64_t> getMatches(QueryMoleculeType &mol, const FingerprintType &fingerprint,
                                     int maxResults, bool &stopFlag) {
        std::vector<uint64_t> results;
        search(mol, fingerprint, maxResults, stopFlag, *_root, results);
        return results;
    }

    std::vector<uint64_t> getMatches(QueryMoleculeType &mol, int maxResults, bool &stopFlag) {
        FingerprintType queryFingerprint(mol);
        return getMatches(mol, queryFingerprint, maxResults, stopFlag);
    }

    explicit Trie(
            std::unique_ptr<std::vector<std::pair<std::unique_ptr<MoleculeType>, std::unique_ptr<FingerprintType>>>> &&data,
            size_t maxDepth = std::numeric_limits<size_t>::max())
            : _maxDepth(maxDepth) {
        _root = std::make_unique<Node>(false);
        buildTrie(*_root, *data);
    }
};
