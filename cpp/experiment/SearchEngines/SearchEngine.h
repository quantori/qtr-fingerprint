#pragma once

#include <vector>
#include <string>


class SearchEngine {
    virtual std::vector<uint64_t> search(const std::string& querySmiles, int maxResults, bool& stopFlag) = 0;
};

