#pragma once

#include <string>

class SearchQuery {
public:
    SearchQuery(SearchQuery &) = default;
    SearchQuery(const SearchQuery&) = default;

    explicit SearchQuery(std::string smiles, size_t maxResults, bool& stopFlag);

    [[nodiscard]] const std::string &smiles() const;

    [[nodiscard]] size_t maxResults() const;

    [[nodiscard]] bool checkStopFlag() const;

    [[nodiscard]] const bool& stopFlag() const;

    void stop();

protected:
    std::string _smiles;
    size_t _maxResults;
    bool& _stopFlag;
};

