#pragma once

#include <string>

#include "answer_filtering/PropertiesFilter.h"
#include "BallTreeDriveSearchEngine.h"
#include "modes/web/IdConverter.h"
#include "BaseLibrary.h"

namespace qtr {

    class SearchData {
    public:
        virtual ~SearchData() = default;

        SearchData(size_t ansCount, size_t threadCount, double timeLimit, bool verificationStage);

        struct Query {
            std::unique_ptr<std::string> smiles;
            std::unique_ptr<Fingerprint> fingerprint;
            BaseLibrary baseLibrary = BaseLibrary::BadOption;

            [[nodiscard]] Fingerprint getFingerprint() const;
        };

        virtual std::unique_ptr<QueryData<CIDType>>
        search(const SearchData::Query &query, const PropertiesFilter::Bounds &queryBounds) = 0;

        virtual BaseLibrary getBaseLibrary() const;

    protected:
        size_t ansCount;
        size_t threadsCount;
        double timeLimit;
        bool verificationStage;
    };

} // qtr
