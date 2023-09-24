#include "BallTreeQueryData.h"

#include <algorithm>

using namespace std;

namespace qtr {

    BallTreeQueryData::BallTreeQueryData(size_t stopAnswersCount, double timeLimit, const IndigoFingerprint &query,
                                         unique_ptr <ByIdAnswerFilter> &&filter) :
            QueryData<CIDType>(stopAnswersCount, timeLimit, std::move(filter)),
            _queryFingerprint(query) {}

    const IndigoFingerprint &BallTreeQueryData::getQueryFingerprint() const {
        return _queryFingerprint;
    }

} // namespace qtr
