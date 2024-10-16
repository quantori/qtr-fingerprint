#include <random>
#include <algorithm>

#include "BallTreeRAMSearchEngine.h"
#include "fingerprint_table_io/FingerprintTableReader.h"
#include "Profiling.h"

namespace qtr {

    std::map<std::string, std::vector<std::pair<size_t, std::filesystem::path>>>
    BallTreeRAMSearchEngine::groupedByDriveLeafFiles() const {
        std::map<std::string, std::vector<std::pair<size_t, std::filesystem::path>>> result;
        for (size_t leafId: getLeafIds()) {
            auto leafFilePath = getLeafDir(leafId);
            result[leafFilePath.parent_path().string()].emplace_back(leafId, leafFilePath);
        }
        return result;
    }

    void
    BallTreeRAMSearchEngine::loadLeafFiles(const std::vector<std::pair<size_t, std::filesystem::path>> &leafsList) {
        LOG(INFO) << "Start loading leaf content from drive " << leafsList.front().second.parent_path().parent_path();
        for (size_t i = 0; i < leafsList.size(); i++) {
            auto &[leafId, leafFilePath] = leafsList[i];
            auto &bucket = _buckets[leafNumberById(leafId)];
            FingerprintTableReader reader(leafFilePath / ("data" + fingerprintTableExtension));
            for (const auto &[id, fingerprint]: reader) {
                assert(id < _totalFingerprints);
                bucket.emplace_back(id, fingerprint);
            }
            if (i % 1'00'000 == 0) {
                LOG(INFO) << "Loaded " << i << "/" << leafsList.size()
                          << " ball tree leafs from drive " << leafsList.front().second.parent_path().parent_path();
            }
        }
        LOG(INFO) << "Finish loading leaf content from drive " << leafsList.front().second.parent_path().parent_path();
    }

    std::vector<CIDType> BallTreeRAMSearchEngine::searchInLeaf(size_t leafId, const Fingerprint &query) const {
        ProfileScope("searchInLeaf");
        std::vector<CIDType> answers;
        for (const auto &[id, fingerprint]: _buckets[leafNumberById((leafId))]) {
            if (query <= fingerprint)
                answers.emplace_back(id);
        }
        return answers;
    }

    std::vector<std::vector<uint64_t>>
    BallTreeRAMSearchEngine::divideLeavesIntoGroups(const std::vector<uint64_t> &leaves, size_t threads) const {
        ProfileScope("divideLeavesIntoGroups");
        auto shuffledLeaves = leaves;
        std::mt19937 random_generator(0);
        std::shuffle(shuffledLeaves.begin(), shuffledLeaves.end(), random_generator);
        std::vector<std::vector<uint64_t>> result(threads);
        for (size_t i = 0; i < shuffledLeaves.size(); i++) {
            result[i % result.size()].emplace_back(shuffledLeaves[i]);
        }
        return std::move(result);
    }

    std::vector<fingerprint_table_value_t> BallTreeRAMSearchEngine::getLeafContent(size_t leafId) const {
        return _buckets.at(leafNumberById((leafId)));
    }

} // qtr
