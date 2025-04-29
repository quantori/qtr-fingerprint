#pragma once

#include <concepts>
#include <vector>

#include "frameworks/FrameworkInterface.h"
#include "dataset/SmilesStorage.h"

template<typename DatasetT>
concept DatasetInterface = requires(DatasetT dataset) {
    typename DatasetT::FrameworkT;
    requires FrameworkInterface<typename DatasetT::FrameworkT>;
    requires std::constructible_from<DatasetT, SmilesStorage &&>;
    {
    dataset.size()
    } -> std::same_as<size_t>;
    {
    dataset.molecule(std::declval<size_t>())
    } -> std::same_as<std::unique_ptr<typename DatasetT::FrameworkT::MoleculeT>>;
};
