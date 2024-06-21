#pragma once

#include <vector>
#include <filesystem>
#include <future>
#include <atomic>
#include <type_traits>


#include "BallTree.h"
#include "Fingerprint.h"
#include "BallTreeTypes.h"
#include "query_data/QueryDataWithFingerprint.h"
#include "answer_filtering/AnswerFilter.h"
#include "answer_filtering/AlwaysTrueFilter.h"
#include "fingerprint_table_io/FingerprintTableIOConsts.h"

namespace qtr {

    class BallTreeSearchEngine : public BallTree {
    public:
        template<typename BinaryReader>
        BallTreeSearchEngine(BinaryReader &nodesReader,
                             std::vector<std::filesystem::path> dataDirectories,
                             size_t fingerprintLength, uint64_t fingerprintsCount);

        void search(QueryDataWithFingerprint &queryData, size_t threads) const;

        virtual ~BallTreeSearchEngine() = default;

        [[nodiscard]] virtual std::vector<fingerprint_table_value_t> getLeafContent(size_t leafId) const = 0;

        [[nodiscard]] std::vector<size_t> getLeafIds() const;

    protected:
        [[nodiscard]] const std::filesystem::path &getLeafDir(size_t nodeId) const;

        void initLeafDataPaths();

        template<typename BinaryReader>
        void loadNodes(BinaryReader &reader);

        void findLeaves(const Fingerprint &fingerprint, size_t currentNode, std::vector<CIDType> &leaves) const;

        [[nodiscard]] virtual std::vector<CIDType>
        searchInLeaf(size_t leafId, const Fingerprint &query) const = 0;

        [[nodiscard]] virtual std::vector<std::vector<uint64_t>>
        divideLeavesIntoGroups(const std::vector<uint64_t> &leaves, size_t threads) const = 0;

        void processLeafGroup(QueryDataWithFingerprint &queryData, const std::vector<uint64_t> &leaves) const;

        void
        oneThreadSearch(QueryDataWithFingerprint &queryData, size_t currentNode,
                        ByIdAnswerFilter &filterObject, bool &stopFlag) const;

        void processLeaf(QueryDataWithFingerprint &queryData, uint64_t leafId, ByIdAnswerFilter &filterObject) const;

        std::vector<std::filesystem::path> _leafDirPaths;
        uint64_t _totalFingerprints;
    };

    template<typename BinaryReader>
    BallTreeSearchEngine::BallTreeSearchEngine(BinaryReader &nodesReader,
                                               std::vector<std::filesystem::path> dataDirectories,
                                               size_t fingerprintLength, uint64_t fingerprintsCount)
            : BallTree(dataDirectories, fingerprintLength), _totalFingerprints(fingerprintsCount) {
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
            _nodes.back().centroid.load(reader, _fingerprintLength);
        }
    }

} // qtr

