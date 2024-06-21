#pragma once

#include "SearchData.h"

namespace qtr {

    class BallTreeSearchData : public SearchData {
    public:
        BallTreeSearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree,
                           std::shared_ptr<const IdConverter> idConverter, size_t ansCount, size_t threadCount,
                           double timeLimit, bool verificationStage);

        std::unique_ptr<QueryData<CIDType>>
        search(const SearchData::Query &query, const PropertiesFilter::Bounds &queryBounds) override;

        std::shared_ptr<const BallTreeSearchEngine> ballTree;
        std::shared_ptr<const IdConverter> idConverter;

    private:
        [[nodiscard]] std::unique_ptr<ByIdAnswerFilter>
        getFilter(const SearchData::Query &query, const PropertiesFilter::Bounds &queryBounds) const;
    };
} // qtr
