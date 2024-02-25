#pragma once

#include <string>

#include "Fingerprint.h"

namespace qtr {

    using fingerprint_table_value_t = typename std::pair<uint64_t, Fingerprint>;

    const std::string fingerprintTableExtension = ".ft";

} // qtr