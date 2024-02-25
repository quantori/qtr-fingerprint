#include <unordered_map>
#include <string>

#include "BallTreeDriveSearchEngine.h"
#include "fingerprint_table_io/FingerprintTableReader.h"

using namespace std;

namespace qtr {
    vector <CIDType> BallTreeDriveSearchEngine::searchInLeaf(size_t leafId, const Fingerprint &query) const {
        vector<CIDType> result;
        for (const auto &[id, fingerprint]: FingerprintTableReader(getFingerprintTablePath(leafId))) {
            if (query <= fingerprint)
                result.emplace_back(id);
        }
        return result;
    }

    filesystem::path BallTreeDriveSearchEngine::getFingerprintTablePath(size_t nodeId) const {
        return getLeafDir(nodeId) / ("data" + fingerprintTableExtension);
    }

    vector <vector<uint64_t>>
    BallTreeDriveSearchEngine::divideLeavesIntoGroups(const vector <uint64_t> &leaves, size_t threads) const {
        unordered_map<string, size_t> groupId;
        vector<vector<uint64_t>> result;
        for (auto &leaf: leaves) {
            auto leafLocation = getLeafDir(leaf).parent_path().string();
            auto it = groupId.find(leafLocation);
            if (it == groupId.end()) {
                it = groupId.emplace(leafLocation, groupId.size()).first;
                result.emplace_back();
            }
            result[it->second].emplace_back(leaf);

        }
        return result;
    }
} // qtr