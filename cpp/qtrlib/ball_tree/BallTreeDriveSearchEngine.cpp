#include "BallTreeDriveSearchEngine.h"
#include "fingerprint_table_io/FingerprintTableReader.h"

namespace qtr {
    void BallTreeDriveSearchEngine::searchInLeaf(size_t leafId, BallTreeSearchEngine::QueryData &queryData) const {
        FingerprintTableReader reader(getLeafFile(leafId));
        for (const auto &[id, fingerprint]: reader) {
            if (queryData.query <= fingerprint)
                putAnswer(id, queryData);
        }
    }
} // qtr