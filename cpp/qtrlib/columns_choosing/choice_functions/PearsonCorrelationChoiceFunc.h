#pragma once

#include "ColumnsChooser.h"

namespace qtr {

    class PearsonCorrelationChoiceFunc {
    public:
        choice_result_t operator()(choice_argumnent_t fingerprints) const;
    };


} // namespace qtr
