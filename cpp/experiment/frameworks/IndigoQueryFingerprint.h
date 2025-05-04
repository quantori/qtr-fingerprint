#pragma once

#include <vector>
#include "base_cpp/array.h"

template<typename T>
class IndigoQueryFingerprint {
public:
    explicit IndigoQueryFingerprint(const indigo::Array<T> &fingerprint) {
        for (int i = 0; i < fingerprint.size(); i++) {
            byte b = fingerprint.at(i);
            if (b == 0) {
                continue;
            }
            _items.emplace_back(i, b);
        }
    }

    [[nodiscard]] bool isSubFingerprint(const indigo::Array<T> &fingerprint) const {
        for (auto &[idx, b]: _items) {
            byte val = fingerprint.at(idx);
            if ((b & val) != b) {
                return false;
            }
        }
        return true;
    }

private:
    std::vector<std::pair<int, T>> _items;
};

