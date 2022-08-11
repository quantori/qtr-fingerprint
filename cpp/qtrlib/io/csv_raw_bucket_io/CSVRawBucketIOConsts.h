#pragma once

#include <utility>
#include <string>

#include "Fingerprint.h"

namespace qtr {

    using csv_raw_bucket_value_t = typename std::pair<std::string, qtr::IndigoFingerprint>;

    const std::string csvFileExtension = ".csv";

    const char csvSplitSymbol = '~';

} // namespace qtr