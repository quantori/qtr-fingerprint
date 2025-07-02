#pragma once

#include "StatTable.h"

struct NodeStat {
    size_t visits = 0;
    size_t skips = 0;

    [[nodiscard]] inline StatRow toStatRow() const {
        StatRow row;
        row.addEntry("visits", visits);
        row.addEntry("skips", skips);
        return row;
    };
};

class BallTreeNodesStat {
private:
    std::vector<NodeStat> _nodeStats;

public:
    inline explicit BallTreeNodesStat(size_t nodesCount) : _nodeStats(nodesCount) {}

    inline const NodeStat &operator[](size_t nodeId) const {
        return _nodeStats[nodeId];
    }

    inline NodeStat &operator[](size_t nodeId) {
        return _nodeStats[nodeId];
    }

    [[nodiscard]] inline StatTable toStatTable() const {
        StatTable table;
        for (size_t i = 0; i < _nodeStats.size(); i++) {
            auto row = _nodeStats[i].toStatRow();
            row.addEntry("nodeId", i);
            table.addRow(row);
        }
        return table;
    }
};