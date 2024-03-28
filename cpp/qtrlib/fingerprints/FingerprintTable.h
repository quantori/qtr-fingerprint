#pragma once

#include "Fingerprint.h"

#include <vector>

namespace qtr {

    class FingerprintTable : public std::vector<Fingerprint> {
    public:
        FingerprintTable() = default;

        explicit FingerprintTable(const std::string &sdfFile);
    };

} // namespace qtr
