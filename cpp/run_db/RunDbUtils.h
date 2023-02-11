#pragma once

#include "HuffmanSmilesTable.h"
#include "BallTreeSearchEngine.h"
#include "PropertiesFilter.h"

namespace qtr {
    std::pair<bool, std::unique_ptr<BallTreeQueryData>>
    doSearch(const std::string &querySmiles, const qtr::BallTreeSearchEngine &ballTree,
             const std::shared_ptr<const SmilesTable> &smilesTable, size_t ansCount, size_t threadsCount,
             const PropertiesFilter::Bounds &queryBounds,
             const std::shared_ptr<const std::vector<PropertiesFilter::Properties>> &propertiesTable);

    std::shared_ptr<SmilesTable>
    loadSmilesTable(const std::filesystem::path &smilesTablePath, const HuffmanCoder &huffmanCoder);


} // namespace qtr