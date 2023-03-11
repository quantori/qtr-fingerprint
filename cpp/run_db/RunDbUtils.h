#pragma once

#include "HuffmanSmilesTable.h"
#include "BallTreeSearchEngine.h"
#include "PropertiesFilter.h"
#include "search_data/SearchData.h"

namespace qtr {

    enum class DbType {
        OnDrive,
        InRam
    };

    std::pair<bool, std::unique_ptr<BallTreeQueryData>>
    runSearch(const SearchData &searchData, const std::string &querySmiles, const PropertiesFilter::Bounds &queryBounds);

} // namespace qtr