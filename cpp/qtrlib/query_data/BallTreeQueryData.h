#pragma once

#include <future>
#include <vector>
#include <memory>

#include "Fingerprint.h"
#include "AlwaysTrueFilter.h"
#include "AnswerFilter.h"
#include "QueryData.h"

namespace qtr {

    class BallTreeQueryData : public QueryData<CIDType> {
    public:
        BallTreeQueryData(size_t stopAnswersCount, double timeLimit,
                          const IndigoFingerprint &query = IndigoFingerprint(),
                          std::unique_ptr<ByIdAnswerFilter> &&filter = std::make_unique<AlwaysTrueFilter<CIDType>>());

        [[nodiscard]] const IndigoFingerprint &getQueryFingerprint() const;

    private:
        IndigoFingerprint _queryFingerprint;
    };

} // namespace qtr
