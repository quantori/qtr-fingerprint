#pragma once

#include "ColumnsSelector.h"

namespace qtr {

    class PrimaryOrderSelectionFunction {
        PrimaryOrderSelectionFunction() = default;

        selection_result_t operator()(selection_argument_t fingerprints) const;
    };

} // namespace qtr
