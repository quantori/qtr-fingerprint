#pragma once

#include <vector>

#include "frameworks/FrameworkInterface.h"

template<typename FrameworkT> requires FrameworkInterface<FrameworkT>
class SearchResult {
public:
    SearchResult() = default;

    void addResult(const FrameworkT::MoleculeT &molecule) {
        _results.emplace_back(molecule);
    }

    [[nodiscard]] const FrameworkT::MoleculeT get(size_t n) const {
        return _results.at(n);
    }

    [[nodiscard]] size_t size() const {
        return _results.size();
    }

    virtual ~SearchResult() = default;

protected:
    std::vector<typename FrameworkT::MoleculeT> _results;
};
