#pragma once

#include <concepts>
#include <cstdlib>

#include "GraphMol/GraphMol.h"

template<typename T, typename Mol>
concept Fingerprint = requires(T t, size_t index) {
    requires std::constructible_from<T, const Mol &>;
    requires std::constructible_from<T>;
    { t.isSubFingerprintOf(t) } -> std::convertible_to<bool>;
    { t.getBit(index) } -> std::convertible_to<bool>;
    { t.size() } -> std::convertible_to<size_t>;
    { t |= t };
};
