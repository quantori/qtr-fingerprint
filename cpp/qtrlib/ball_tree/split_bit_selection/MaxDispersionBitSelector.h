#pragma once

#include "split_bit_selection/BitSelector.h"

namespace qtr {

    class MaxDispersionBitSelector : public BitSelector {
    public:
        size_t operator()(const ColumnsStatistic &statistic) const override;
    };

} // qtr

