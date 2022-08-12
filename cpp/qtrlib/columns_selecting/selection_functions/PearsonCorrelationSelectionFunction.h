#pragma once

#include "ColumnsSelector.h"

namespace qtr {
    class PearsonCorrelationSelectionFunction {
    private:
        static const size_t subsetSize = 1000;

    public:
        selection_result_t operator()(select_argument_t fingerprints) const;
    };


} // namespace qtr
