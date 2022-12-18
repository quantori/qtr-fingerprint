#pragma once

#include <vector>
#include <filesystem>
#include <future>
#include <functional>

#include "BallTree.h"
#include "Fingerprint.h"


namespace qtr {
    using CIDType = uint64_t;

    namespace {
        auto noFiltering = [](CIDType) {
            return true;
        };
    }

    class BallTreeSearchEngine : public BallTree {
    public:
        template<typename BinaryReader>
        BallTreeSearchEngine(BinaryReader &nodesReader, std::vector<std::filesystem::path> dataDirectories);

        struct QueryData {
            const IndigoFingerprint &query;
            std::vector<CIDType> &result;
            std::mutex &resultLock;
            size_t ansCount;
            bool isTerminate;
            const std::function<bool(CIDType)> &filter;

            void updateIsTerminate();

            void addAnswer(CIDType value);
        };

    protected:
        void initLeafDataPaths();

        void searchInSubtree(size_t nodeId, QueryData &queryData) const;

        [[nodiscard]] const std::filesystem::path &getLeafFile(size_t nodeId) const;

        template<typename BinaryReader>
        void loadNodes(BinaryReader &reader);

        [[nodiscard]] std::vector<size_t> getLeafIds() const;

        static void putAnswer(CIDType ansValue, QueryData &queryData);

        void findLeafs(const IndigoFingerprint &fingerprint, size_t currentNode, std::vector<CIDType> &leafs) const;

        virtual void searchInLeaf(size_t leafId, QueryData &queryData) const = 0;

        void processLeafGroup(QueryData &queryData, std::vector<uint64_t> leafs, size_t group,
                              size_t totalGroups) const;
    public:
        std::vector<CIDType> search(const IndigoFingerprint &query, size_t ansCount, size_t startDepth,
                                    const std::function<bool(CIDType)> &filter = noFiltering) const;

    protected:
        std::vector<std::filesystem::path> _leafDataPaths;
    };

    template<typename BinaryReader>
    BallTreeSearchEngine::BallTreeSearchEngine(BinaryReader &nodesReader,
                                               std::vector<std::filesystem::path> dataDirectories)
            : BallTree(dataDirectories) {
        loadNodes(nodesReader);
        assert(__builtin_popcountll(_nodes.size() + 1) == 1);
        _depth = log2Floor(_nodes.size());
        assert((1ull << (_depth + 1)) == _nodes.size() + 1);
        initLeafDataPaths();
    }

    template<typename BinaryReader>
    void BallTreeSearchEngine::loadNodes(BinaryReader &reader) {
        uint64_t treeSize;
        reader.read((char *) &treeSize, sizeof treeSize);
        _nodes.reserve(treeSize);
        for (size_t i = 0; i < treeSize; i++) {
            _nodes.emplace_back();
            _nodes.back().centroid.load(reader);
        }
    }

} // qtr

