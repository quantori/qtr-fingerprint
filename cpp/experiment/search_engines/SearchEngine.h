#pragma once

#include <vector>
#include <string>


class SearchEngine {
public:
    virtual std::vector<uint64_t> getMatches(const std::string &querySmiles, int maxResults, bool &stopFlag) = 0;

    virtual ~SearchEngine() = default;
};

