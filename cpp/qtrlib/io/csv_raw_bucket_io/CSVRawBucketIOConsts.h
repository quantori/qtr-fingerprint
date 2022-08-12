#pragma once

#include <utility>
#include <string>
#include <sstream>

#include "Fingerprint.h"

namespace qtr {

    using csv_raw_bucket_value_t = typename std::pair<std::string, qtr::IndigoFingerprint>;

    const std::string csvFileExtension = ".csv";

    const std::string csvRawBucketIndexColName = "SMILES";

    const char csvSplitSymbol = '~';

    inline std::string csvRawBucketHeader() {
        std::stringstream stream;
        stream << csvRawBucketIndexColName << csvSplitSymbol;
        for (size_t i = 0; i < IndigoFingerprint::sizeInBits; i++) {
            stream << i;
            stream << (i + 1 == IndigoFingerprint::sizeInBits ? '\n' : csvSplitSymbol);
        }
        return stream.str();
    }

    const size_t csvRawBucketHeaderSize = csvRawBucketHeader().size();

} // namespace qtr