#include "SearchResult.h"

void SearchResult::addResultByIndex(size_t resIdx) {
    _resultIndices.push_back(resIdx);
}

size_t SearchResult::getNthResult(size_t n) const {
    return _resultIndices.at(n);
}

size_t SearchResult::size() const {
    return _resultIndices.size();
}
