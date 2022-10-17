#pragma once

#include <vector>
#include <filesystem>
#include <future>
#include <functional>

#include "BallTree.h"
#include "Fingerprint.h"


namespace qtr {

    namespace {
        auto noFiltering = [](size_t) {
            return true;
        };
    }

    class BallTreeSearchEngine : public BallTree {
    public:
        template<typename BinaryReader>
        BallTreeSearchEngine(BinaryReader &nodesReader, std::vector<std::filesystem::path> dataDirectories);

    protected:
        void initLeafDataPaths();

        struct QueryData {
            const IndigoFingerprint &query;
            std::vector<size_t> &result;
            std::mutex &resultLock;
            size_t ansCount;
            bool isTerminate;
            const std::function<bool(size_t)> &filter;

            void updateIsTerminate();

            void addAnswer(size_t value);
        };

        void searchInSubtree(size_t nodeId, QueryData &queryData) const;

        const std::filesystem::path &getLeafFile(size_t nodeId) const;

        std::filesystem::path &getLeafFile(size_t nodeId);

        template<typename BinaryReader>
        void loadNodes(BinaryReader &reader);

        std::vector<size_t> getLeafIds() const;

        static void putAnswer(size_t ansValue, QueryData& queryData);

        virtual void searchInLeaf(size_t leafId, QueryData &queryData) const = 0;

    public:
        std::vector<size_t> search(const IndigoFingerprint &query, size_t ansCount, size_t startDepth,
                                   const std::function<bool(size_t)> &filter = noFiltering) const;

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

