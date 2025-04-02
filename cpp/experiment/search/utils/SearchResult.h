#pragma once

#include <vector>


class SearchResult {
public:
    SearchResult() = default;

    void addResultByIndex(size_t resIdx);

    [[nodiscard]] size_t getNthResult(size_t n) const;

    [[nodiscard]] size_t size() const;

    virtual ~SearchResult() = default;

protected:
    std::vector<size_t> _resultIndices;
};
