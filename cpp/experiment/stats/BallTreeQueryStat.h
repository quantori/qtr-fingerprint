#pragma once

#include "StatRow.h"
#include "Profiling.h"

struct BallTreeQueryStat {
    size_t leafSearches = 0;
    std::vector<size_t> nodesVisitedPerDepth;
    std::vector<size_t> subsetSizePerDepth;

    inline explicit BallTreeQueryStat(size_t treeDepth) : nodesVisitedPerDepth(treeDepth + 1, 0),
                                                          subsetSizePerDepth(treeDepth + 1, 0) {
    }

    [[nodiscard]] inline StatRow toStatRow() const {
        ProfileScope("BallTreeQueryStat::toStatRaw");
        StatRow row;
        row.addEntry("leafSearches", leafSearches);
        for (size_t i = 0; i < nodesVisitedPerDepth.size(); i++) {
            row.addEntry("nodesVisitedAtDepth" + std::to_string(i), nodesVisitedPerDepth[i]);
        }
        for (size_t i = 0; i < subsetSizePerDepth.size(); i++) {
            row.addEntry("subsetSizeAtDepth" + std::to_string(i), subsetSizePerDepth[i]);
        }
        return row;
    }
};