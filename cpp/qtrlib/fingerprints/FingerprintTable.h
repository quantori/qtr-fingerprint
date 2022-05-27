#pragma once

#include "Fingerprint.h"

#include <cstdint>

/**
 * @tparam FINGERPRINT_SIZE fingerprint size in bits
 */
template<std::size_t FINGERPRINT_SIZE>
class FingerprintTable {
public:
    /**
     * Create new table of fingerprints with given size
     * @param size
     */
    FingerprintTable(std::size_t size) {
        _data = new Fingerprint<FINGERPRINT_SIZE>[size];
    }

    ~FingerprintTable() {
        delete _data;
    }

    Fingerprint<FINGERPRINT_SIZE> &operator[](std::size_t index) {
        return _data[index];
    }

    const Fingerprint<FINGERPRINT_SIZE> &operator[](std::size_t index) const {
        return _data[index];
    }

private:
    Fingerprint<FINGERPRINT_SIZE> *_data = nullptr;
};
