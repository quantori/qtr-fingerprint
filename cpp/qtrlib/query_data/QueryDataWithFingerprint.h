#pragma once

#include <future>
#include <vector>
#include <memory>

#include "Fingerprint.h"
#include "AlwaysTrueFilter.h"
#include "AnswerFilter.h"
#include "QueryData.h"

namespace qtr {

    class QueryDataWithFingerprint : public QueryData<CIDType> {
    public:
        QueryDataWithFingerprint(size_t stopAnswersCount, double timeLimit,
                                 const Fingerprint &query,
                                 std::unique_ptr<ByIdAnswerFilter> &&filter = std::make_unique<AlwaysTrueFilter<CIDType>>());

        [[nodiscard]] const Fingerprint &getQueryFingerprint() const;

    private:
        Fingerprint _queryFingerprint;
    };

} // namespace qtr
