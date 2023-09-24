#pragma once

#include "SearchData.h"

namespace qtr {

    class QtrSearchData : public SearchData {
    public:
        QtrSearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree,
                      std::shared_ptr<const IdConverter> idConverter,
                      size_t ansCount, size_t threadCount, double timeLimit);

        std::unique_ptr<QueryData<CIDType>>
        search(const std::string &querySmiles, const PropertiesFilter::Bounds &queryBounds) override;


        std::shared_ptr<const BallTreeSearchEngine> ballTree;
        std::shared_ptr<const IdConverter> idConverter;
    };

} // qtr
