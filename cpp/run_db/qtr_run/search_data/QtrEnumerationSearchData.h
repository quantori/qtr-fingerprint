#pragma once

#include "SearchData.h"
#include "CFStorage.h"

namespace qtr {

    class QtrEnumerationSearchData : public SearchData {
    public:
        QtrEnumerationSearchData(std::shared_ptr<std::vector<Fingerprint>> allFingerprints, size_t ansCount,
                                 size_t threadCount, double timeLimit,
                                 bool verificationStage);

        std::unique_ptr<QueryData<CIDType>>
        search(const SearchData::Query &query, const PropertiesFilter::Bounds &) override;

    private:
        std::shared_ptr<std::vector<Fingerprint>> _allFingerprints;
    };

} // qtr

