#pragma once

#include <concepts>
#include <cstdlib>

#include "GraphMol/GraphMol.h"

template<typename T>
concept Fingerprint = requires(T t, size_t index) {
    requires std::constructible_from<T, const typename T::MoleculeType &>;
    requires std::constructible_from<T>;
    { t.isSubFingerprintOf(t) } -> std::convertible_to<bool>;
    { t.getBit(index) } -> std::convertible_to<bool>;
    { t.size() } -> std::convertible_to<size_t>;
    { t |= t };
};

//template<typename T, typename BaseFingerprintT>
//concept QueryFingerprint = requires(T t, const BaseFingerprintT& fp) {
//    requires std::constructible_from<T, const BaseFingerprintT&>;
//    { t.isSubFingerprintOf(fp) } -> std::convertible_to<bool>;
//};
