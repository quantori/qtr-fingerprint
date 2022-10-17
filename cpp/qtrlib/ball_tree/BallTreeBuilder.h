#pragma once

#include <vector>
#include <filesystem>
#include <mutex>

#include "Fingerprint.h"
#include "BallTree.h"
#include "split_bit_selection/BitSelector.h"

namespace qtr {

    class BallTreeBuilder : public BallTree {
    private:

        void buildTree(size_t depth, size_t parallelizationDepth, const BitSelector &bitSelector);

        std::vector<size_t> buildFirstLevels(size_t nodeId, size_t levels, const BitSelector &bitSelector);

        void buildLastLevels(size_t nodeId, size_t levels, const BitSelector &bitSelector,
                             const std::filesystem::path &dataDirectory);

        void splitNodeManyFiles(size_t nodeId, const BitSelector &bitSelector);

        void
        splitNodeOneFile(size_t nodeId, const BitSelector &bitSelector, const std::filesystem::path &dataDirectory);

        std::vector<std::filesystem::path> getNodeFiles(size_t nodeId) const;

        void calculateCentroid(size_t nodeId);

        void deleteNodeFromFilesystem(size_t nodeId) const;

    public:
        BallTreeBuilder(size_t depth, size_t parallelizationDepth, std::vector<std::filesystem::path> dataDirectories,
                        const BitSelector &bitSelector);

        template<typename BinaryWriter>
        void dumpNodes(BinaryWriter &writer);
    };

    template<typename BinaryWriter>
    void BallTreeBuilder::dumpNodes(BinaryWriter &writer) {
        uint64_t treeSize = _nodes.size();
        writer.write((char *) &treeSize, sizeof treeSize);
        for (auto &_node: _nodes) {
            _node.centroid.dump(writer);
        }
    }

} // qtr
