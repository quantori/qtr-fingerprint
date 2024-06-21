#pragma once

#include <concepts>
#include <cstdlib>

template<typename T>
concept Fingerprint = requires(T t, size_t index) {
    { t.isSubFingerprintOf(t) } -> std::convertible_to<bool>;
    { t.getBit(index) } -> std::convertible_to<bool>;
};
