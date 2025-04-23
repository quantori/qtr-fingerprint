#pragma once

#include <vector>

#include "frameworks/FrameworkInterface.h"

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

    virtual ~SearchResult() = default;

protected:
    std::vector<ResultT> _results;
};
