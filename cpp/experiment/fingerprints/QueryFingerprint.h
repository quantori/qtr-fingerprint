#pragma once

#include <vector>

#include "FingerprintConcept.h"

template<Fingerprint FP>
class QueryFingerprint {
public:
    explicit QueryFingerprint(const FP& baseFingerprint): _bitPositions() {
        for (size_t i = 0; i < baseFingerprint.size(); i++) {
            if (baseFingerprint.getBit(i)) {
                _bitPositions.emplace_back(i);
            }
        }
    }

    bool isSubFingerprintOf(const FP& other) const {
        for (size_t i : _bitPositions) {
            if (!other.getBit(i)) {
                return false;
            }
        }
        return true;
    }

    [[nodiscard]] const std::vector<size_t>& bits() const {
        return _bitPositions;
    }

private:
    std::vector<size_t> _bitPositions;
};