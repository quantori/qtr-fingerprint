#pragma once

#include "SmilesTable.h"
#include "BallTreeSearchEngine.h"

namespace qtr {
    std::pair<bool, std::unique_ptr<BallTreeQueryData>> doSearch(const std::string &querySmiles,
                                                                 const qtr::BallTreeSearchEngine &ballTree,
                                                                 const std::shared_ptr<const SmilesTable> &smilesTable,
                                                                 size_t ansCount);

    SmilesTable loadSmilesTable(const std::filesystem::path &smilesTablePath, const HuffmanCoder &huffmanCoder);


} // namespace qtr