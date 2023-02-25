#include "BallTreeDriveSearchEngine.h"
#include "fingerprint_table_io/FingerprintTableReader.h"

using namespace std;

namespace qtr {
    vector <CIDType> BallTreeDriveSearchEngine::searchInLeaf(size_t leafId, const IndigoFingerprint &query) const {

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
} // qtr