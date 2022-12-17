#include <utility>
#include <vector>
#include "BallTree.h"


namespace qtr {

    BallTree::BallTree(std::vector<std::filesystem::path> dataDirectories)
            : _nodes(), _depth(0), _dataDirectories(std::move(dataDirectories)) {}


    size_t BallTree::leftChild(size_t nodeId) {
        return nodeId * 2 + 1;
    }

    size_t BallTree::rightChild(size_t nodeId) {
        return nodeId * 2 + 2;
    }

    size_t BallTree::parent(size_t nodeId) {
        return (nodeId - 1) >> 1ull;
    }

    bool BallTree::isLeaf(size_t nodeId) const {
        return leftChild(nodeId) >= _nodes.size();
    }

    size_t BallTree::root() {
        return 0;
    }

    size_t BallTree::leafNumberById(size_t leafId) const {
        assert((1ull << _depth) - 1 <= leafId);
        leafId -= (1ull << _depth) - 1;
        assert(leafId < (1ull << _depth));
        return leafId;
    }

} // qtr