#pragma once

#include <string>

#include "Fingerprint.h"

namespace qtr {

    using fingerprint_table_value_t = typename std::pair<size_t, IndigoFingerprint>;

    const std::string fingerprintTableExtension = ".ft";

} // qtr