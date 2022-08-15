#include <cassert>
#include <random>
#include <utility>

#include "PearsonCorrelationSelectionFunction.h"
#include "PearsonCorrelationSelectionFunctionUtils.h"
#include "ColumnsSelectionUtils.h"
#include "FingerprintTable.h"


namespace qtr {

    selection_result_t
    PearsonCorrelationSelectionFunction::operator()(const IndigoFingerprintTable &fingerprints) const {
        assert(!fingerprints.empty());
        IndigoFingerprintTable fingerprintSubset;
        if (_fingerprintSubsetSize != 0) {
            fingerprintSubset = selectFingerprintSubset(fingerprints, _fingerprintSubsetSize);
        } else {
            fingerprintSubset = fingerprints;
        }
        std::map<size_t, std::vector<bool>> columns = fingerprintsToColumns(fingerprintSubset, _columnsSubset);
        std::map<size_t, double> maxCorrelation = findMaxAbsPearsonCorrelation(columns);
        std::vector<size_t> columnsIndexes = sortIndexesByValues(maxCorrelation);
        return columnsIndexes;
    }

    PearsonCorrelationSelectionFunction::PearsonCorrelationSelectionFunction(std::vector<size_t> columnsSubset,
                                                                             size_t fingerprintSubsetSize)
            : _columnsSubset(std::move(columnsSubset)), _fingerprintSubsetSize(fingerprintSubsetSize) {}

} // namespace qtr
