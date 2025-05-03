#include "SearchQuery.h"

#include <utility>

const std::string &SearchQuery::smiles() const {
    return _smiles;
}

SearchQuery::SearchQuery(std::string smiles, size_t maxResults, bool &stopFlag) : _smiles(std::move(smiles)),
                                                                                  _maxResults(maxResults),
                                                                                  _stopFlag(stopFlag) {
}

bool SearchQuery::checkStopFlag() const {
    return _stopFlag;
}

void SearchQuery::stop() {
    _stopFlag = true;
}

size_t SearchQuery::maxResults() const {
    return _maxResults;
}

const bool &SearchQuery::stopFlag() const {
    return _stopFlag;
}

