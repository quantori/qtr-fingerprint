#include "QueryDataWithFingerprint.h"

#include <algorithm>

using namespace std;

namespace qtr {

    QueryDataWithFingerprint::QueryDataWithFingerprint(size_t stopAnswersCount, double timeLimit, const Fingerprint &query,
                                                       unique_ptr <ByIdAnswerFilter> &&filter) :
            QueryData<CIDType>(stopAnswersCount, timeLimit, std::move(filter)),
            _queryFingerprint(query) {}

    const Fingerprint &QueryDataWithFingerprint::getQueryFingerprint() const {
        return _queryFingerprint;
    }

} // namespace qtr
