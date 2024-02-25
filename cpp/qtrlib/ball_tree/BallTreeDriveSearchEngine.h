#pragma once

#include "BallTreeSearchEngine.h"

namespace qtr {

    class BallTreeDriveSearchEngine : public BallTreeSearchEngine {
    protected:
        [[nodiscard]] std::filesystem::path getFingerprintTablePath(size_t nodeId) const;

        [[nodiscard]] std::vector<std::vector<uint64_t>>
        divideLeavesIntoGroups(const std::vector<uint64_t> &leaves, size_t threads) const override;

    public:

        using BallTreeSearchEngine::BallTreeSearchEngine;

        [[nodiscard]] std::vector<CIDType> searchInLeaf(size_t leafId, const Fingerprint &query) const override;
    };

} // qtr