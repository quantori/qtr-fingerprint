#pragma once

#include "BallTreeSearchEngine.h"

#include <map>

#include "fingerprint_table_io/FingerprintTableReader.h"

namespace qtr {

    class BallTreeRAMSearchEngine : public BallTreeSearchEngine {
    public:
        template<typename BinaryReader>
        BallTreeRAMSearchEngine(BinaryReader &nodesReader,
                                std::vector<std::filesystem::path> dataDirectories,
                                size_t fingerprintLength, size_t totalFingerprints);

        [[nodiscard]] std::vector<CIDType> searchInLeaf(size_t leafId, const Fingerprint &query) const override;

        [[nodiscard]] std::vector<fingerprint_table_value_t> getLeafContent(size_t leafId) const override;


    protected:
        [[nodiscard]] std::map<std::string, std::vector<std::pair<size_t, std::filesystem::path>>>
        groupedByDriveLeafFiles() const;

        void loadLeafFiles(const std::vector<std::pair<size_t, std::filesystem::path>> &leafsList);

        [[nodiscard]] std::vector<std::vector<uint64_t>>
        divideLeavesIntoGroups(const std::vector<uint64_t> &leaves, size_t threads) const override;

        std::vector<std::vector<fingerprint_table_value_t>> _buckets;
    };

    template<typename BinaryReader>
    BallTreeRAMSearchEngine::BallTreeRAMSearchEngine(BinaryReader &nodesReader,
                                                     std::vector<std::filesystem::path> dataDirectories,
                                                     size_t fingerprintLength, size_t totalFingerprints)
            :BallTreeSearchEngine(nodesReader, dataDirectories, fingerprintLength, totalFingerprints),
             _buckets(1ull << _depth) {
        std::vector<std::future<void>> tasks;
        auto groupedLeafFiles = groupedByDriveLeafFiles();
        for (const auto &[_, leafs]: groupedLeafFiles) {
            tasks.emplace_back(
                    std::async(std::launch::async, &BallTreeRAMSearchEngine::loadLeafFiles, this, std::cref(leafs))
            );
        }
        for (auto &task: tasks) {
            task.get();
        }
    }

} // qtr
