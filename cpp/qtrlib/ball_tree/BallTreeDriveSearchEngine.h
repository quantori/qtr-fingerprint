#pragma once

#include "BallTreeSearchEngine.h"

namespace qtr {

    class BallTreeDriveSearchEngine : public BallTreeSearchEngine {
    protected:
        [[nodiscard]] std::filesystem::path getFingerprintTablePath(size_t nodeId) const;

    public:
        [[nodiscard]] std::vector<CIDType> searchInLeaf(size_t leafId, const IndigoFingerprint &query) const override;
    };

} // qtr