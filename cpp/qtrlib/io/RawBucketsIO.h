#pragma once

#include <utility>
#include <string>

#include "Fingerprint.h"

namespace qtr {

    using bucket_value_t = typename std::pair<std::string, qtr::IndigoFingerprint>;

} // namespace qtr

#include "RawBucketReader.h"
#include "RawBucketWriter.h"