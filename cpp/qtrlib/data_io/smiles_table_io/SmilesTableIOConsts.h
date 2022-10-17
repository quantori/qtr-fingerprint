#pragma once

#include <utility>
#include <string>

namespace qtr {
    using smiles_table_value_t = std::pair<uint64_t, std::string>;

    const std::string smilesTableExtension = ".st";
}