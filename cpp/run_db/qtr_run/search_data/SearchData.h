#pragma once

#include <string>

#include "answer_filtering/PropertiesFilter.h"
#include "BallTreeDriveSearchEngine.h"
#include "modes/web/IdConverter.h"

namespace qtr {

    class SearchData {
    public:
        virtual ~SearchData() = default;

        SearchData(size_t ansCount, size_t threadCount, double timeLimit, bool verificationStage);

        struct Query {
            std::unique_ptr<std::string> smiles;
            std::unique_ptr<Fingerprint> fingerprint;
        };

        virtual std::unique_ptr<QueryData<CIDType>>
        search(const SearchData::Query &query, const PropertiesFilter::Bounds &queryBounds) = 0;

    protected:
        size_t ansCount;
        size_t threadsCount;
        double timeLimit;
        bool verificationStage;
    };

} // qtr
