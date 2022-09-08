#pragma once

#include <vector>
#include <filesystem>

#include "Fingerprint.h"
#include "split_bit_selectors/BitSelector.h"

namespace qtr {

    class BallTree {
    public:
        class Node;

    private:
        std::vector<Node> _nodes;
        std::vector<std::filesystem::path> _dataDirectories;
        std::vector<std::filesystem::path> _leafDataPath;

        static size_t leftChild(size_t nodeId);

        static size_t rightChild(size_t nodeId);

        static size_t parent(size_t nodeId);

        bool isLeaf(size_t nodeId) const;

        size_t root() const;

        void build();

        std::vector<size_t> _buildFirstLevels(size_t nodeId, size_t levels, const BitSelector& bitSelector);

        void _splitNodeManyFiles(size_t nodeId, size_t splitBit);

        std::vector<std::filesystem::path> _getNodesFiles(size_t nodeId) const;

        void prepareFilesForNode(size_t nodeId) const;

    public:
        BallTree(size_t depth, size_t parallelizationDepth, std::vector<std::filesystem::path> dataDirectories,
                 const BitSelector &bitSelector);

        BallTree() = default;

        template<typename BinaryWriter>
        void dump(BinaryWriter &writer) {
            // todo
        }

        template<typename BinaryReader>
        void load(BinaryReader &reader) {
            // todo
        }
    };

    class BallTree::Node {
    private:
        IndigoFingerprint centroid;
    };

} // qtr
