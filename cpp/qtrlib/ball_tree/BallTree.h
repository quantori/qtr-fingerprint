#pragma once

#include <vector>
#include <filesystem>

#include "Fingerprint.h"

namespace qtr {

    class BallTree {
    public:
        class Node;

    protected:
        explicit BallTree(std::vector<std::filesystem::path> dataDirectories, size_t fingerprintLength);

        static size_t leftChild(size_t nodeId);

        static size_t rightChild(size_t nodeId);

        static size_t parent(size_t nodeId);

        [[nodiscard]] bool isLeaf(size_t nodeId) const;

        [[nodiscard]] size_t leafNumberById(size_t leafId) const;

        static size_t root();

    protected:
        std::vector<Node> _nodes;
        size_t _depth;
        std::vector<std::filesystem::path> _dataDirectories;
        size_t _fingerprintLength;
    };

    struct BallTree::Node {
        Fingerprint centroid;
    };

} // qtr
