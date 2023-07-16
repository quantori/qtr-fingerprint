#pragma once

#include <vector>
#include <filesystem>
#include <future>
#include <atomic>
#include <type_traits>


#include "BallTree.h"
#include "Fingerprint.h"
#include "BallTreeTypes.h"
#include "query_data/BallTreeQueryData.h"
#include "answer_filtering/AnswerFilter.h"
#include "answer_filtering/AlwaysTrueFilter.h"

namespace qtr {

    class BallTreeSearchEngine : public BallTree {
    public:

        inline static std::atomic<double> ballTreeSearchTimer = 0;

        template<typename BinaryReader>
        BallTreeSearchEngine(BinaryReader &nodesReader, std::vector<std::filesystem::path> dataDirectories);

        void search(BallTreeQueryData &queryData, size_t threads) const;

        virtual ~BallTreeSearchEngine() = default;

    protected:
        [[nodiscard]] const std::filesystem::path &getLeafDir(size_t nodeId) const;

        void initLeafDataPaths();

        template<typename BinaryReader>
        void loadNodes(BinaryReader &reader);

        [[nodiscard]] std::vector<size_t> getLeafIds() const;

        void findLeaves(const IndigoFingerprint &fingerprint, size_t currentNode, std::vector<CIDType> &leaves) const;

        [[nodiscard]] virtual std::vector<CIDType> searchInLeaf(size_t leafId, const IndigoFingerprint &query) const = 0;

        [[nodiscard]] virtual std::vector<std::vector<uint64_t>>
        divideLeavesIntoGroups(const std::vector<uint64_t> &leaves, size_t threads) const = 0;

        void processLeafGroup(BallTreeQueryData &queryData, const std::vector <uint64_t> &leaves) const;

        std::vector<std::filesystem::path> _leafDirPaths;
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

