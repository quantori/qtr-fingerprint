#pragma once

#include "ColumnsSelector.h"

namespace qtr {

    class IotaOrderSelectionFunction {
    public:
        selection_result_t operator()(selection_argument_t fingerprints) const;

        IotaOrderSelectionFunction() = default;
    };

} // namespace qtr
