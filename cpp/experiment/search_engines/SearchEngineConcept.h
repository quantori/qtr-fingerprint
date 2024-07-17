#pragma once

#include <vector>
#include <concepts>

#include "FingerprintConcept.h"

enum class SearchEngineType {
    QtrRDKit,
    QtrIndigo,
    RDKit,
    Indigo,
};

template<typename T>
concept SearchEngine = requires(T t, const std::string &querySmiles, int maxResults, bool &stopFlag,
                                const T::QueryMoleculeType &queryMolecule, const T::FingerprintType& queryFingerprint,
                                const T::StorageMoleculeType& storageMolecule, const T::MoleculeType& molecule) {
    { t.getMatches(queryMolecule, maxResults, stopFlag) } -> std::convertible_to<std::vector<uint64_t>>;
    { t.getMatches(queryMolecule, queryFingerprint, maxResults, stopFlag) } -> std::convertible_to<std::vector<uint64_t>>;

    { t.smilesToMolecule(querySmiles) } -> std::convertible_to<std::unique_ptr<typename T::MoleculeType>>;
    { t.smilesToQueryMolecule(querySmiles) } -> std::convertible_to<std::unique_ptr<typename T::QueryMoleculeType>>;
    { t.storageMoleculeToMolecule(storageMolecule) } -> std::convertible_to<std::unique_ptr<typename T::MoleculeType>>;
    { t.moleculeToStorageMolecule(molecule) } -> std::convertible_to<std::unique_ptr<typename T::StorageMoleculeType>>;

    requires std::constructible_from<T, const std::vector<std::string> &>;
    requires std::constructible_from<T, std::vector<std::unique_ptr<typename T::StorageMoleculeType>> &&>;
    requires std::constructible_from<T, std::vector<std::pair<std::unique_ptr<typename T::StorageMoleculeType>, std::unique_ptr<typename T::FingerprintType>>> &&>;
};

