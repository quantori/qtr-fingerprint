#pragma once

#include <utility>
#include <string>

namespace qtr {
    using string_table_value_t = std::pair<uint64_t, std::string>;

    const std::string stringTableExtension = ".st";
}