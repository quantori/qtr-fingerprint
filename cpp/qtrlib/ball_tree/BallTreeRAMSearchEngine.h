#pragma once

#include "BallTreeSearchEngine.h"

#include <map>

#include "fingerprint_table_io/FingerprintTableReader.h"

#include "IndigoQueryMolecule.h"
#include "IndigoSubstructureMatcher.h"

namespace qtr {

    class BallTreeRAMSearchEngine : public BallTreeSearchEngine {
    public:
        template<typename BinaryReader>
        BallTreeRAMSearchEngine(BinaryReader &nodesReader,
                                std::vector<std::filesystem::path> dataDirectories);

        virtual void searchInLeaf(size_t leafId, QueryData &queryData) const override;
        virtual void searchInLeafIds(const indigo_cpp::IndigoSessionPtr &indigoSessionPtr,
                                     const indigo_cpp::IndigoQueryMolecule &queryMol,
                                     size_t leafId, QueryData &queryData) const override;
    protected:
        std::map<std::string, std::vector<std::pair<size_t, std::filesystem::path>>> groupedByDriveLeafFiles() const;

        void loadLeafFiles(const std::vector<std::pair<size_t, std::filesystem::path>> &leafsList);

    protected:
        std::vector<std::vector<fingerprint_table_value_t>> _buckets;
    };

    template<typename BinaryReader>
    BallTreeRAMSearchEngine::BallTreeRAMSearchEngine(BinaryReader &nodesReader,
                                                     std::vector<std::filesystem::path> dataDirectories)
            :BallTreeSearchEngine(nodesReader, dataDirectories), _buckets(1ull << _depth) {
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
