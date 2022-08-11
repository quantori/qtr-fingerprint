#include <cassert>
#include <numeric>
#include <random>

#include "PearsonCorrelationChoiceFunc.h"
#include "PearsonCorrelationChoiceFuncUtils.h"
#include "ColumnsChoiceUtils.h"
#include "FingerprintTable.h"


namespace qtr {

    choice_result_t PearsonCorrelationChoiceFunc::operator()(const IndigoFingerprintTable &fingerprints) const {
        assert(!fingerprints.empty());
        auto subset = chooseSubset(fingerprints, subsetSize);
        auto columns = fingerprintsToColumns(subset);
        auto maxCorrelation = findMaxAbsPearsonCorrelation(columns);
        auto columnsIndexes = sortIndexesByValues(maxCorrelation);
        return columnsIndexes;
    }

} // namespace qtr
