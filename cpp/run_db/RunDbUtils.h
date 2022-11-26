#pragma once

#include "SmilesTable.h"
#include "BallTreeDriveSearchEngine.h"

namespace qtr {
    std::pair<bool, std::vector<uint64_t>>
    doSearch(const std::string &querySmiles, const BallTreeSearchEngine &ballTree,
                    const SmilesTable &smilesTable, uint64_t ansCount, uint64_t startSearchDepth);

    HuffmanCoder buildHuffmanCoder(const std::filesystem::path &smilesTablePath);

    SmilesTable loadSmilesTable(const std::filesystem::path &smilesTablePath, const HuffmanCoder &huffmanCoder);

    std::pair<HuffmanCoder, SmilesTable> loadCoderAndTable(const std::filesystem::path &smilesTablePath);

} // namespace qtr