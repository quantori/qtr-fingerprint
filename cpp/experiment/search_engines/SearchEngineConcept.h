#pragma once

#include <vector>
#include <concepts>

#include "FingerprintConcept.h"

enum class SearchEngineType {
    BallTreeRDKit,
    BallTreeIndigo,
    TrieRDKit,
    TrieIndigo,
    RDKit,
    Indigo,
    IndigoBruteForce,
};

template<typename T>
concept SearchEngine = requires(T t, const std::string &querySmiles, int maxResults, bool &stopFlag,
        T::QueryMoleculeType &queryMolecule, T::FingerprintType& queryFingerprint,
        T::StorageMoleculeType& storageMolecule, T::MoleculeType& molecule) {
    { t.getMatches(queryMolecule, maxResults, stopFlag) } -> std::convertible_to<std::vector<uint64_t>>;
    { t.getMatches(queryMolecule, queryFingerprint, maxResults, stopFlag) } -> std::convertible_to<std::vector<uint64_t>>;

    { t.smilesToMolecule(querySmiles) } -> std::convertible_to<std::unique_ptr<typename T::MoleculeType>>;
    { t.smilesToQueryMolecule(querySmiles) } -> std::convertible_to<std::unique_ptr<typename T::QueryMoleculeType>>;
    { t.storageMoleculeToMolecule(storageMolecule) } -> std::convertible_to<std::unique_ptr<typename T::MoleculeType>>;
    { t.moleculeToStorageMolecule(molecule) } -> std::convertible_to<std::unique_ptr<typename T::StorageMoleculeType>>;

    requires std::constructible_from<T, std::unique_ptr<std::vector<std::string>> &&>;
    requires std::constructible_from<T, std::unique_ptr<std::vector<std::unique_ptr<typename T::StorageMoleculeType>>> &&>;
    requires std::constructible_from<T, std::unique_ptr<std::vector<std::pair<std::unique_ptr<typename T::StorageMoleculeType>, std::unique_ptr<typename T::FingerprintType>>>> &&>;
};

