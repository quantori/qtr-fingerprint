#pragma once

#include "SearchData.h"

namespace qtr {

    class RDKitSearchData : public SearchData {
    public:
        RDKitSearchData(const std::filesystem::path &dbDataDir, size_t ansCount, size_t threadsCount,
                        double timeLimit, bool verificationStage);

        std::unique_ptr<QueryData<CIDType>>
        search(const SearchData::Query &query, const PropertiesFilter::Bounds &) override;

        ~RDKitSearchData() override;

    private:
        std::shared_ptr<RDKit::SubstructLibrary> _substructLibrary;
    };

} // qtr
