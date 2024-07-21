#pragma once

#include <glog/logging.h>

#include "FingerprintConcept.h"
#include "SearchEngineConcept.h"

#include "GraphMol/SmilesParse/SmilesParse.h"

#include <string>

template<SearchEngine SE>
class BallTree {
//    static_assert(std::is_same<typename FP::MoleculeType, typename SE::MoleculeType>::value);
// TODO: think about static asserts
    using MoleculeType = typename SE::MoleculeType;
    using QueryMoleculeType = typename SE::QueryMoleculeType;
    using StorageMoleculeType = typename SE::StorageMoleculeType;
    using FingerprintType = typename SE::FingerprintType;
private:
    struct Node {
        FingerprintType centroid;
    };

    std::vector<Node> _nodes;
    std::vector<std::unique_ptr<SE>> _leafSearchEngines;
    size_t _depth;

    static size_t leftChild(size_t nodeId) {
        return nodeId * 2 + 1;
    }

    static size_t rightChild(size_t nodeId) {
        return nodeId * 2 + 2;
    }

    static size_t parent(size_t nodeId) {
        return (nodeId - 1) >> 1;
    }

    [[nodiscard]] bool isLeaf(size_t nodeId) const {
        return leftChild(nodeId) >= _nodes.size();
    }

    [[nodiscard]] size_t nodeIdToLeafId(size_t leafId) const {
        assert((1ull << _depth) - 1 <= leafId);
        leafId -= (1ull << _depth) - 1;
        assert(leafId < (1ull << _depth));
        return leafId;
    }

    static size_t root() {
        return 0;
    }

    static size_t endNodeId() {
        return -1;
    }

    [[nodiscard]] size_t traverseUp(size_t nodeId) const {
        while (true) { // TODO: could be optimized with bit operations
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

    void
    searchInLeafNode(const QueryMoleculeType &mol, const FingerprintType &fingerprint, int maxResults, bool &stopFlag,
                     size_t nodeId,
                     std::vector<uint64_t> &results) {
        assert(isLeaf(nodeId));
        size_t leafId = nodeIdToLeafId(nodeId);
        std::vector<uint64_t> leafResults = _leafSearchEngines[leafId]->getMatches(mol, fingerprint,
                                                                                   maxResults - results.size(),
                                                                                   stopFlag);
        results.insert(results.end(), leafResults.begin(), leafResults.end());
    }

    static std::vector<size_t>
    countBitStat(
            const std::vector<std::pair<std::unique_ptr<StorageMoleculeType>, std::unique_ptr<FingerprintType>>> &data) {
        std::vector<size_t> bitStat(data[0].second->size());
        for (auto &[mol, fp]: data) {
            assert(fp->size() == bitStat.size());
            for (size_t i = 0; i < bitStat.size(); i++) {
                bitStat[i] += fp->getBit(i);
            }
        }
        return bitStat;
    }

    static size_t getBalance(const std::vector<size_t> &bitStat, size_t i) {
        return std::abs((int) bitStat.size() - 2 * (int) bitStat[i]);
    }

    static size_t selectSplitBit(const std::vector<size_t> &bitStat) {
        assert(!bitStat.empty());
        size_t resBit = 0;
        for (size_t i = 0; i < bitStat.size(); i++) {
            if (getBalance(bitStat, resBit) > getBalance(bitStat, i)) {
                resBit = i;
            }
        }
        return resBit;
    }

    static void splitData(
            std::unique_ptr<std::vector<std::pair<std::unique_ptr<StorageMoleculeType>, std::unique_ptr<FingerprintType>>>> &&data,
            std::vector<std::pair<std::unique_ptr<StorageMoleculeType>, std::unique_ptr<FingerprintType>>> &leftData,
            std::vector<std::pair<std::unique_ptr<StorageMoleculeType>, std::unique_ptr<FingerprintType>>> &rightData) {
        auto bitStat = countBitStat(*data);
        auto splitBit = selectSplitBit(bitStat);
        for (auto &[mol, fp]: *data) {
            auto &destData = leftData.size() * 2 >= data->size()
                             || (rightData.size() * 2 < data->size() && fp->getBit(splitBit) == 1) ? rightData
                                                                                                   : leftData;
            destData.emplace_back(std::move(mol), std::move(fp));
        }
        assert(leftData.size() * 2 <= data->size() + 1);
        assert(rightData.size() * 2 <= data->size() + 1);
    }


    static FingerprintType evaluateCentroid(
            const std::vector<std::pair<std::unique_ptr<StorageMoleculeType>, std::unique_ptr<FingerprintType >>> &data) {
        FingerprintType res;
        for (auto &[mol, fp]: data) {
            res |= *fp;
        }
        return res;
    }

public:
    std::vector<uint64_t>
    getMatches(const QueryMoleculeType &mol, const FingerprintType &fingerprint, int maxResults, bool &stopFlag) {
        size_t nodeId = root();
        std::vector<uint64_t> results;
        while (nodeId != endNodeId()
               && !stopFlag
               && (int) results.size() < maxResults) {
            auto &node = _nodes[nodeId];
            bool shouldSkipSubtree = !fingerprint.isSubFingerprintOf(node.centroid);
            if (!shouldSkipSubtree && isLeaf(nodeId)) {
                searchInLeafNode(mol, fingerprint, maxResults, stopFlag, nodeId, results);
            }
            if (shouldSkipSubtree) {
                nodeId = traverseUp(nodeId);
            } else {
                nodeId = traverseDown(nodeId);
            }
        }
        return results;
    }

    std::vector<uint64_t> getMatches(const QueryMoleculeType &mol, int maxResults, bool &stopFlag) {
        FingerprintType queryFingerprint(mol);
        return getMatches(mol, queryFingerprint, maxResults, stopFlag);
    }

    explicit BallTree(
            std::unique_ptr<std::vector<std::pair<std::unique_ptr<StorageMoleculeType>, std::unique_ptr<FingerprintType>>>> &&data) {
        using NodeDataType = std::vector<std::pair<std::unique_ptr<StorageMoleculeType>, std::unique_ptr<FingerprintType>>>;
//        _depth = 1;
        _depth = std::max((size_t) 2, (size_t) std::ceil(std::log2(data->size() / 25))); // TODO: 25 - magic constant
        _leafSearchEngines.resize(1ull << _depth); // TODO: check +- 1
        _nodes.resize((2ull << _depth) - 1);
        std::vector<std::unique_ptr<NodeDataType>> nodesData((1u << _depth) * 2 - 1);
        nodesData[root()] = std::move(data);
        for (size_t nodeId = root(); nodeId < nodesData.size(); nodeId++) {
            if (isLeaf(nodeId)) {
                size_t leafId = nodeIdToLeafId(nodeId);
                _nodes[nodeId].centroid = evaluateCentroid(*nodesData[nodeId]);
                _leafSearchEngines[leafId] = std::make_unique<SE>(std::move(nodesData[nodeId]));
            } else {
                nodesData[leftChild(nodeId)] = std::make_unique<NodeDataType>();
                nodesData[rightChild(nodeId)] = std::make_unique<NodeDataType>();
                splitData(std::move(nodesData[nodeId]), *nodesData[leftChild(nodeId)], *nodesData[rightChild(nodeId)]);
            }
        }
        for (int nodeId = (1u << _depth) - 2; nodeId >= 0; nodeId--) {
            assert(!isLeaf(nodeId));
            _nodes[nodeId].centroid |= _nodes[leftChild(nodeId)].centroid;
            _nodes[nodeId].centroid |= _nodes[rightChild(nodeId)].centroid;
        }
//        LOG(INFO) << "Ball Tree:";
//        for (size_t i = 0; i < _nodes.size(); i++) {
//            if (!isLeaf(i)) {
//                assert(_nodes[leftChild(i)].centroid.isSubFingerprintOf(_nodes[i].centroid));
//                assert(_nodes[rightChild(i)].centroid.isSubFingerprintOf(_nodes[i].centroid));
//            }
//            std::stringstream ss;
//            ss << i << ": ";
//            for (size_t j = 0; j < _nodes[i].centroid.size(); j++) {
//                ss << _nodes[i].centroid.getBit(j);
//            }
//            LOG(INFO) << ss.str();
//        }
    }
};