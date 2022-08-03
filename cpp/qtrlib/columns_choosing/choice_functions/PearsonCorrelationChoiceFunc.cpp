#include <cassert>

#include "PearsonCorrelationChoiceFunc.h"
#include "PearsonCorrelationChoiceFuncUtils.h"
#include "ColumnsChoiceUtils.h"
#include "FingerprintTable.h"


namespace qtr {

    choice_result_t PearsonCorrelationChoiceFunc::operator()(const IndigoFingerprintTable &fingerprints) const {
        assert(!fingerprints.empty());
        auto columns = fingerprintsToColumns(fingerprints);
        auto maxCorrelation = findMaxAbsPearsonCorrelation(columns);
        auto columnsIndexes = sortIndexesByValues(maxCorrelation);
        return columnsIndexes;
    }

} // namespace qtr
