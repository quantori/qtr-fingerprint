#pragma once

#include <cstdlib>
#include <filesystem>
#include <vector>


namespace qtr {

    class BitSelector {
    public:
        virtual size_t operator()(const std::vector<std::filesystem::path> &fpTablePaths) const = 0;

        virtual size_t operator()(const std::filesystem::path &fpTablePath) const = 0;
    };

} // qtr

