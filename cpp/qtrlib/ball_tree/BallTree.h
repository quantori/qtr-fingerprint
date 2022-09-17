#pragma once

#include <vector>
#include <filesystem>
#include <mutex>

#include "Fingerprint.h"
#include "split_bit_selection//BitSelector.h"

namespace qtr {

    class BallTree {
    public:
        class Node;

    private:
        std::vector<Node> _nodes;
        std::vector<std::filesystem::path> _dataDirectories;
        std::vector<std::filesystem::path> _leafDataPaths;
        size_t _depth;

        static size_t leftChild(size_t nodeId);

        static size_t rightChild(size_t nodeId);

        static size_t parent(size_t nodeId);

        bool isLeaf(size_t nodeId) const;

        size_t root() const;

        void buildTree(size_t depth, size_t parallelizationDepth, const BitSelector &bitSelector);

        std::vector<size_t> buildFirstLevels(size_t nodeId, size_t levels, const BitSelector &bitSelector);

        void buildLastLevels(size_t nodeId, size_t levels, const BitSelector &bitSelector,
                             const std::filesystem::path &dataDirectory);

        void splitNodeManyFiles(size_t nodeId, const BitSelector &bitSelector);

        void
        splitNodeOneFile(size_t nodeId, const BitSelector &bitSelector, const std::filesystem::path &dataDirectory);

        std::vector<std::filesystem::path> getNodeFiles(size_t nodeId) const;

        void calculateCentroid(size_t nodeId);

        void initLeafDataPaths();

        void
        searchInSubtree(size_t nodeId, const IndigoFingerprint &query, size_t ansCount, std::vector<size_t> &result,
                        std::mutex &resultLock, bool &isTerminate) const;

        const std::filesystem::path &getLeafFile(size_t nodeId) const;

    public:
        BallTree(size_t depth, size_t parallelizationDepth, std::vector<std::filesystem::path> dataDirectories,
                 const BitSelector &bitSelector);

        template<typename BinaryReader>
        BallTree(BinaryReader &nodesReader, std::vector<std::filesystem::path> dataDirectories);

        template<typename BinaryWriter>
        void dumpNodes(BinaryWriter &writer);

        template<typename BinaryReader>
        void loadNodes(BinaryReader &reader);

        std::vector<size_t> search(const IndigoFingerprint &query, size_t ansCount, size_t startDepth) const;
    };

    struct BallTree::Node {
        IndigoFingerprint centroid;
    };

    template<typename BinaryReader>
    BallTree::BallTree(BinaryReader &nodesReader, std::vector<std::filesystem::path> dataDirectories)
            : _dataDirectories(std::move(dataDirectories)) {
        loadNodes(nodesReader);
        assert(__builtin_popcountll(_nodes.size() + 1) == 1);
        _depth = log2Floor(_nodes.size());
        assert((1ull << (_depth + 1)) == _nodes.size() + 1);
        initLeafDataPaths();
    }

    template<typename BinaryWriter>
    void BallTree::dumpNodes(BinaryWriter &writer) {
        uint64_t treeSize = _nodes.size();
        writer.write((char *) &treeSize, sizeof treeSize);
        for (auto &_node: _nodes) {
            _node.centroid.dump(writer);
        }
    }

    template<typename BinaryReader>
    void BallTree::loadNodes(BinaryReader &reader) {
        uint64_t treeSize;
        reader.read((char *) &treeSize, sizeof treeSize);
        _nodes.reserve(treeSize);
        for (size_t i = 0; i < treeSize; i++) {
            _nodes.emplace_back();
            _nodes.back().centroid.load(reader);
        }
    }

} // qtr
