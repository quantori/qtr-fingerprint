#pragma once

#include "split_bit_selectors/BitSelector.h"

namespace qtr {

    class MaxDispersionBitSelector : public BitSelector {
    public:
        size_t operator()(const std::vector<std::filesystem::path> &fpTablePaths) const override;

        size_t operator()(const std::filesystem::path &fpTablePath) const override;
    };

} // qtr

