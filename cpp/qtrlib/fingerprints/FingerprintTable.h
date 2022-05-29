#pragma once

#include "Fingerprint.h"

#include <cstdint>

class FingerprintTableTestFixture;

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
    FingerprintTable(std::size_t size) : _size(size) {
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


    struct Iterator {
        Iterator(Fingerprint<FINGERPRINT_SIZE> *data) : _data(data) {}

        Fingerprint<FINGERPRINT_SIZE> &operator*() const {
            return *_data;
        }

        Fingerprint<FINGERPRINT_SIZE> *operator->() {
            return _data;
        }

        Iterator& operator++() {
            ++_data;
            return *this;
        }

        Iterator operator++(int) {
            Iterator tmp = _data;
            ++_data;
            return tmp;
        }

        friend bool operator==(const Iterator& a, const Iterator& b) {
            return a._data == b._data;
        }

        friend bool operator!=(const Iterator& a, const Iterator& b) {
            return a._data != b._data;
        }
    private:
        Fingerprint<FINGERPRINT_SIZE> *_data;
    };

    Iterator begin() {
        return Iterator(_data);
    }

    Iterator end() {
        return Iterator(_data + _size);
    }

private:
    Fingerprint<FINGERPRINT_SIZE> *_data = nullptr;
    const std::size_t _size;

    friend FingerprintTableTestFixture;
};
