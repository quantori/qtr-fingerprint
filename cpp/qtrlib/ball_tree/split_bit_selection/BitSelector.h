#pragma once

#include <cstdlib>
#include <filesystem>
#include <vector>

#include "ball_tree/split_bit_selection/ColumnsStatistic.h"


namespace qtr {

    class BitSelector {
    public:
        virtual size_t operator()(const ColumnsStatistic &statistic) const = 0;

        virtual ~BitSelector() = default;
    };

} // qtr

