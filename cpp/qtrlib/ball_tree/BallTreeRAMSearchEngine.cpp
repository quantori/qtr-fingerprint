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
        LOG(INFO) << "Start loading leaf content from drive " << leafsList.front().second.parent_path().parent_path();
        for (auto &[leafId, leafFilePath]: leafsList) {
            FingerprintTableReader reader(leafFilePath);
            std::copy(reader.begin(), reader.end(), std::back_inserter(_buckets[leafNumberById(leafId)]));
        }
        LOG(INFO) << "Finish loading leaf content from drive " << leafsList.front().second.parent_path().parent_path();
    }

    std::vector<CIDType> BallTreeRAMSearchEngine::searchInLeaf(size_t leafId, const IndigoFingerprint &query) const {
        std::vector<CIDType> answers;
        for (const auto &[id, fingerprint]: _buckets[leafNumberById((leafId))]) {
            if (query <= fingerprint)
                answers.emplace_back(id);
        }
        return answers;
    }

} // qtr
