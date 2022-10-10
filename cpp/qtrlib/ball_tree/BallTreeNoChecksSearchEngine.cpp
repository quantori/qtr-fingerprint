#include "BallTreeNoChecksSearchEngine.h"
#include "fingerprint_table_io/FingerprintTableReader.h"

namespace qtr {

    void BallTreeNoChecksSearchEngine::loadLeafData() {
        auto leafIds = getLeafIds();
        _leafData.resize(leafIds.size());
        for (auto &leafId: getLeafIds()) {
            std::filesystem::path leafFile = _leafDataPaths[leafNumberById(leafId)];
            auto &leafData = _leafData[leafNumberById(leafId)];
            for (const auto &[fingerprintId, _]: FingerprintTableReader(leafFile)) {
                leafData.emplace_back(fingerprintId);
            }
        }
    }

    void BallTreeNoChecksSearchEngine::searchInLeaf(size_t leafId, BallTreeSearchEngine::QueryData &queryData) const {
        for (auto& fingerprintId : _leafData[leafNumberById(leafId)]) {
            putAnswer(fingerprintId, queryData);
        }
    }

} // qtr