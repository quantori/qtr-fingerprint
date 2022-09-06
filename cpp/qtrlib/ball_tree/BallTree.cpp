#include <utility>
#include <future>

#include "BallTree.h"

namespace qtr {

    namespace {

        void splitFileBybit(const std::filesystem::path &file, size_t splitBit, const std::filesystem::path &leftFile,
                            const std::filesystem::path &rightFiles) {
            // todo как сплить ровно по ровну? как синхронизироваться между файлами?
        }
    } // namespace

    BallTree::BallTree(size_t depth, size_t parallelizationDepth, std::vector<std::filesystem::path> dataDirectories,
                       const BitSelector &bitSelector)
            : _dataDirectories(std::move(dataDirectories)) {
        _nodes.resize((1ull << depth) - 1);
        _leafDataPath.resize((1ull << (depth - 1)));
        build();
        // todo
    }

    size_t BallTree::leftChild(size_t nodeId) {
        return nodeId * 2 + 1;
    }

    size_t BallTree::rightChild(size_t nodeId) {
        return nodeId * 2 + 2;
    }

    size_t BallTree::root() const {
        return 0;
    }

    bool BallTree::isLeaf(size_t nodeId) const {
        return leftChild(nodeId) >= _nodes.size();
    }

    void BallTree::build() {
        // todo
    }

    void BallTree::_splitNodeManyFiles(size_t nodeId, size_t splitBit) {
        // todo
        prepareFilesForNode(leftChild(nodeId));
        prepareFilesForNode(rightChild(nodeId));
        std::vector<std::future<void>> tasks;
        tasks.reserve(_dataDirectories.size());
        for (auto &nodeFile: _getNodesFiles(nodeId)) {

        }

    }

    std::vector<size_t> BallTree::_buildFirstLevels(size_t nodeId, size_t levels, const BitSelector &bitSelector) {
        if (levels == 0)
            return {nodeId};
        size_t splitBit = bitSelector(_getNodesFiles(nodeId));
        _splitNodeManyFiles(nodeId, splitBit);
        std::vector<size_t> nodes = _buildFirstLevels(leftChild(nodeId), levels - 1, bitSelector);
        std::vector<size_t> rightNodes = _buildFirstLevels(rightChild(nodeId), levels - 1, bitSelector);
        nodes.insert(nodes.end(), rightNodes.begin(), rightNodes.end());
        return nodes;
    }

    std::vector<std::filesystem::path> BallTree::_getNodesFiles(size_t nodeId) const {
        std::vector<std::filesystem::path> result;
        for (auto &dirPath: _dataDirectories) {
            for (auto &filePath: findFiles(dirPath / std::to_string(nodeId), "")) {
                result.emplace_back(filePath);
            }
        }
        return result;
    }

    size_t BallTree::parent(size_t nodeId) {
        return (nodeId - 1) >> 1ull;
    }

    void BallTree::prepareFilesForNode(size_t nodeId) const {
        for (auto &dirPath: _dataDirectories) {
            std::filesystem::create_directory(dirPath / std::to_string(nodeId));
        }
    }

} // qtr