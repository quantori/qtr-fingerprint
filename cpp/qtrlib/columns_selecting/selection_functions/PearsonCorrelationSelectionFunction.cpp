#include <cassert>
#include <random>

#include "PearsonCorrelationSelectionFunction.h"
#include "PearsonCorrelationSelectionFunctionUtils.h"
#include "ColumnsSelectionUtils.h"
#include "FingerprintTable.h"


namespace qtr {

    selection_result_t
    PearsonCorrelationSelectionFunction::operator()(const IndigoFingerprintTable &fingerprints) const {
        assert(!fingerprints.empty());
//        auto subset = chooseSubset(fingerprints, subsetSize);
        auto subset = fingerprints;
        auto columns = fingerprintsToColumns(subset);
        auto maxCorrelation = findMaxAbsPearsonCorrelation(columns);
        auto columnsIndexes = sortIndexesByValues(maxCorrelation);
        return columnsIndexes;
    }

} // namespace qtr
