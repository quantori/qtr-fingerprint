#pragma once

#include <vector>

#include "FrameworkInterface.h"


template<typename FrameworkT> requires FrameworkInterface<FrameworkT>
class BallTreeFingerprintChecker {
public:
    explicit BallTreeFingerprintChecker(const FrameworkT::FingerprintT &fingerprint) {
        for (size_t bit = 0; bit < FrameworkT::getFingerprintSize(); bit++) {
            if (FrameworkT::getFingerprintBit(fingerprint, bit)) {
                _bits.emplace_back(bit);
            }
        }
    }

    bool check(const FrameworkT::FingerprintT &fingerprint) const {
        for (size_t bit: _bits) {
            if (!FrameworkT::getFingerprintBit(fingerprint, bit)) {
                return false;
            }
        }
        return true;
    }

private:
    std::vector<size_t> _bits;
};