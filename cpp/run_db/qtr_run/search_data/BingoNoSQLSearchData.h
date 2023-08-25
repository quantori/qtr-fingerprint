#pragma once

#include "SearchData.h"
#include "IndigoSession.h"
#include "BingoNoSQL.h"

namespace qtr {

    class BingoNoSQLSearchData : public SearchData {
    public:
        BingoNoSQLSearchData(const std::filesystem::path &dbDataDir, size_t ansCount, size_t threadsCount,
                             double timeLimit);

        std::unique_ptr<QueryData<CIDType>>
        search(const std::string &querySmiles, const PropertiesFilter::Bounds &) override;

        ~BingoNoSQLSearchData() override;

        indigo_cpp::BingoMolecule db;
    };

} // qtr
