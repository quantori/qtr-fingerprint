#pragma once

#include "SmilesTable.h"
#include "BallTreeDriveSearchEngine.h"

namespace qtr {
    std::pair<bool, std::vector<uint64_t>>
    doSearch(const std::string &querySmiles, const BallTreeSearchEngine &ballTree,
                    const SmilesTable &smilesTable, uint64_t ansCount, uint64_t startSearchDepth);

    SmilesTable loadSmilesTable(const std::filesystem::path &smilesTablePath, const HuffmanCoder &huffmanCoder);


} // namespace qtr