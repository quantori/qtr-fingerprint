#include "IotaOrderSelectionFunction.h"

#include <numeric>

namespace qtr {

    selection_result_t IotaOrderSelectionFunction::operator()(const IndigoFingerprintTable &fingerprints) const {
        selection_result_t result(IndigoFingerprint::size);
        std::iota(result.begin(), result.end(), 0);
        return result;
    }

} // namespace qtr
