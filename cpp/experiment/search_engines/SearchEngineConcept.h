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

template<typename T, typename Mol>
concept SearchEngine = requires(T t, const std::string &querySmiles, int maxResults, bool &stopFlag,
                                const Mol &molecule) {
    { t.getMatches(molecule, maxResults, stopFlag) } -> std::convertible_to<std::vector<uint64_t>>;

// TODO: add it later
//    { t.idToSmiles(std::declval<uint64_t>()) } -> std::convertible_to<std::string>;
//    { t.idToMolecule(std::declval<uint64_t>()) } -> std::convertible_to<Mol>;

    { t.smilesToMolecule(querySmiles) } -> std::convertible_to<std::unique_ptr<Mol>>;
// TODO: add it later
//    { t.moleculeToSmiles(molecule) } -> std::convertible_to<std::string>;

    requires std::constructible_from<T, const std::vector<std::string> &>;
    requires std::constructible_from<T, std::vector<std::unique_ptr<Mol>> &&>;
    requires std::constructible_from<T, std::vector<std::pair<std::unique_ptr<Mol>, std::unique_ptr<typename T::FingerprintType>>> &&>;
};
