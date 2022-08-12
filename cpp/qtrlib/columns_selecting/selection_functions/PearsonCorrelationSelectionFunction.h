#pragma once

#include "ColumnsSelector.h"

namespace qtr {
    class PearsonCorrelationSelectionFunction {
    public:
        PearsonCorrelationSelectionFunction(std::vector<size_t> columnsSubset = {},
                                            size_t fingerprintSubsetSize = 0);

        selection_result_t operator()(selection_argument_t fingerprints) const;

    private:
        std::vector<size_t> _columnsSubset;
        size_t _fingerprintSubsetSize;

    };


} // namespace qtr
