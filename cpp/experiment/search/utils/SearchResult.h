#pragma once

#include <utility>
#include <vector>
#include <sstream>

#include "frameworks/FrameworkInterface.h"
#include "StatTable.h"

template<typename ResultT>
class SearchResult {
public:
    SearchResult() = default;

    void addResult(const ResultT &molecule) {
        _results.emplace_back(molecule);
    }

    [[nodiscard]] const ResultT &get(size_t n) const {
        return _results.at(n);
    }

    [[nodiscard]] size_t size() const {
        return _results.size();
    }

    void setStatRow(StatRow &&statRow) {
        _statRow = statRow;
    }

    [[nodiscard]] const StatRow &getStatRow() const {
        return _statRow;
    }

private:
    std::vector<ResultT> _results;
    StatRow _statRow;
};
