#pragma once

#include "SearchData.h"
#include "GraphMol/SubstructLibrary/SubstructLibrary.h"


namespace qtr {

    class RDKitSearchData : public SearchData {
    public:
        RDKitSearchData(const std::filesystem::path &moleculesDir, const std::filesystem::path &fingerprintTablesDir,
                        size_t ansCount, size_t threadsCount, double timeLimit, bool verificationStage);

        std::unique_ptr<QueryData<CIDType>>
        search(const SearchData::Query &query, const PropertiesFilter::Bounds &) override;

        ~RDKitSearchData() override = default;

    private:
        std::shared_ptr<RDKit::SubstructLibrary> _substructLibrary;
    };

} // qtr
