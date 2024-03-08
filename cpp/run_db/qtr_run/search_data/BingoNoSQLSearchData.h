#pragma once

#include "SearchData.h"
#include "IndigoSession.h"
#include "BingoNoSQL.h"

namespace qtr {

    class BingoNoSQLSearchData : public SearchData {
    public:
        BingoNoSQLSearchData(const std::filesystem::path &dbDataDir, size_t ansCount, size_t threadsCount,
                             double timeLimit, bool verificationStage);

        std::unique_ptr<QueryData<CIDType>>
        search(const SearchData::Query &query, const PropertiesFilter::Bounds &) override;

        ~BingoNoSQLSearchData() override;

        indigo_cpp::BingoMolecule db;
    };

} // qtr
