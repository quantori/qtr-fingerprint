#pragma once

#include "SmilesTable.h"
#include "BallTreeDriveSearchEngine.h"

namespace qtr {
    std::pair<bool, std::vector <std::future<void>>>
    doSearch(const std::string &querySmiles, BallTreeDriveSearchEngine::QueryData &queryData,
             const qtr::BallTreeSearchEngine &ballTree,
             const SmilesTable &smilesTable, uint64_t startSearchDepth);

    SmilesTable loadSmilesTable(const std::filesystem::path &smilesTablePath, const HuffmanCoder &huffmanCoder);


} // namespace qtr