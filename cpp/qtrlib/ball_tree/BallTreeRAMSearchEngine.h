#pragma once

#include "BallTreeSearchEngine.h"

#include <map>

#include "fingerprint_table_io/FingerprintTableReader.h"

namespace qtr {

    class BallTreeRAMSearchEngine : public BallTreeSearchEngine {
    public:
        template<typename BinaryReader>
        BallTreeRAMSearchEngine(BinaryReader &nodesReader,
                                std::vector<std::filesystem::path> dataDirectories);

        virtual void searchInLeaf(size_t leafId, QueryData &queryData) const override;

    protected:
        std::map<std::string, std::vector<std::pair<size_t, std::filesystem::path>>> groupedByDriveLeafFiles() const;

        void loadLeafFiles(const std::vector<std::pair<size_t, std::filesystem::path>> &leafsList);

        std::vector<fingerprint_table_value_t> &getLeafBucket(size_t leafId);

        const std::vector<fingerprint_table_value_t> &getLeafBucket(size_t leafId) const;

    protected:
        std::vector<std::vector<fingerprint_table_value_t>> _buckets;
    };

    template<typename BinaryReader>
    BallTreeRAMSearchEngine::BallTreeRAMSearchEngine(BinaryReader &nodesReader,
                                                     std::vector<std::filesystem::path> dataDirectories)
            :BallTreeSearchEngine(nodesReader, dataDirectories), _buckets(1ull << _depth) {
        std::vector<std::future<void>> tasks;
        auto grupedLeafFiles = groupedByDriveLeafFiles();
        for (const auto &[_, leafs]: grupedLeafFiles) {
            tasks.emplace_back(
                    std::async(std::launch::async, &BallTreeRAMSearchEngine::loadLeafFiles, this, std::cref(leafs))
            );
        }
        for (auto &task: tasks) {
            task.get();
        }
    }

} // qtr
