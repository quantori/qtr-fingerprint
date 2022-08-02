#pragma once

#include <utility>
#include <string>

#include "Fingerprint.h"

namespace qtr {

    using raw_bucket_value_t = typename std::pair<std::string, qtr::IndigoFingerprint>;

    const std::string rawBucketExtension = "rb";

} // namespace qtr