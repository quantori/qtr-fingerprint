#include <cassert>
#include <numeric>

#include "PearsonCorrelationChoiceFunc.h"
#include "PearsonCorrelationChoiceFuncUtils.h"
#include "ColumnsChoiceUtils.h"
#include "FingerprintTable.h"


namespace qtr {

    choice_result_t PearsonCorrelationChoiceFunc::operator()(const IndigoFingerprintTable &fingerprints) const {
        assert(!fingerprints.empty());
        std::vector<int> columns(IndigoFingerprint::sizeInBits);
        std::iota(columns.begin(), columns.end(), 0);
        return columns;
//        auto columns = fingerprintsToColumns(fingerprints);
//        auto maxCorrelation = findMaxAbsPearsonCorrelation(columns);
//        auto columnsIndexes = sortIndexesByValues(maxCorrelation);
//        return columnsIndexes;
// todo correlation choice
    }

} // namespace qtr
