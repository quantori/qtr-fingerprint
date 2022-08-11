#pragma once

#include "ColumnsChooser.h"

namespace qtr {
    class PearsonCorrelationChoiceFunc {
    private:
        static const size_t subsetSize = 1000;

    public:
        choice_result_t operator()(choice_argumnent_t fingerprints) const;
    };


} // namespace qtr
