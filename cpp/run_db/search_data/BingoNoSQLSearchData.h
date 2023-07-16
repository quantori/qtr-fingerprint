#pragma once

#include "SearchData.h"
#include "IndigoSession.h"

namespace qtr {

    class BingoNoSQLSearchData : public SearchData {
    public:
        BingoNoSQLSearchData(int db, indigo_cpp::IndigoSessionPtr indigoSessionPtr,
                             TimeTicker &timeTicker, size_t ansCount, size_t threadsCount,
                             double timeLimit);

        std::unique_ptr<QueryData<CIDType>>
        search(const std::string &querySmiles, const PropertiesFilter::Bounds &) override;

        ~BingoNoSQLSearchData() override;


        int db;
        indigo_cpp::IndigoSessionPtr indigoSessionPtr;
    };

} // qtr
