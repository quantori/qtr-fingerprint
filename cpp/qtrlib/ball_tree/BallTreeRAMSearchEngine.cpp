#include "BallTreeRAMSearchEngine.h"

#include "fingerprint_table_io/FingerprintTableReader.h"

namespace qtr {

    std::map<std::string, std::vector<std::pair<size_t, std::filesystem::path>>>
    BallTreeRAMSearchEngine::groupedByDriveLeafFiles() const {
        std::map<std::string, std::vector<std::pair<size_t, std::filesystem::path>>> result;
        for (size_t leafId: getLeafIds()) {
            auto leafFilePath = getLeafFile(leafId);
            result[leafFilePath.parent_path().parent_path().string()].emplace_back(leafId, leafFilePath);
        }
        return result;
    }

    void
    BallTreeRAMSearchEngine::loadLeafFiles(const std::vector<std::pair<size_t, std::filesystem::path>> &leafsList) {
        for (auto &[leafId, leafFilePath]: leafsList) {
            FingerprintTableReader reader(leafFilePath);
            std::copy(reader.begin(), reader.end(), std::back_inserter(getLeafBucket(leafId)));
        }
    }

    const std::vector<fingerprint_table_value_t> &BallTreeRAMSearchEngine::getLeafBucket(size_t leafId) const {
        assert((1ull << _depth) - 1 <= leafId);
        leafId -= (1ull << _depth) - 1;
        assert(leafId <= _buckets.size());
        return _buckets[leafId];
    }

    std::vector<fingerprint_table_value_t> &BallTreeRAMSearchEngine::getLeafBucket(size_t leafId) {
        return const_cast<std::vector<fingerprint_table_value_t> &>(
                const_cast<const BallTreeRAMSearchEngine *>(this)->getLeafBucket(leafId)
        );
    }

    void BallTreeRAMSearchEngine::searchInLeaf(size_t leafId, BallTreeSearchEngine::QueryData &queryData) const {
        for (const auto &[id, fingerprint]: getLeafBucket(leafId)) {
            if (queryData.query <= fingerprint)
                putAnswer(id, queryData);
        }
    }

} // qtr
