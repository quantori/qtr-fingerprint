#pragma once

#include "QtrRamSearchData.h"
#include "CFStorage.h"

namespace qtr {

    class QtrEnumerationSearchData : public QtrRamSearchData {
    public:
        QtrEnumerationSearchData(QtrRamSearchData searchData);

        std::unique_ptr<QueryData<CIDType>>
        search(const SearchData::Query &query, const PropertiesFilter::Bounds &) override;
    };

} // qtr

