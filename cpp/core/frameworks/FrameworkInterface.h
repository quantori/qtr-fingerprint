#pragma once

#include <concepts>
#include <string>
#include <memory>
#include "utils/Config.h"


template<typename FrameworkT>
concept FrameworkInterface = requires(FrameworkT framework) {
    typename FrameworkT::MoleculeT;
    typename FrameworkT::QueryMoleculeT;
    typename FrameworkT::StorageMoleculeT;

    typename FrameworkT::FingerprintT;
    typename FrameworkT::QueryFingerprintT;

//    requires std::is_copy_constructible_v<typename FrameworkT::MoleculeT>;
//    requires std::is_copy_constructible_v<typename FrameworkT::QueryMoleculeT>;
//    requires std::is_copy_constructible_v<typename FrameworkT::StorageMoleculeT>;

    {
    framework.init(std::declval<const Config &>())
    } -> std::same_as<void>;

    {
    FrameworkT::getInstance()
    } -> std::same_as<FrameworkT&>;

    {
    framework.moleculeFromSmiles(std::declval<const std::string &>())
    } -> std::same_as<std::unique_ptr<typename FrameworkT::MoleculeT>>;

    {
    framework.moleculeToSmiles(std::declval<typename FrameworkT::MoleculeT &>())
    } -> std::convertible_to<std::string>;

    {
    framework.queryMoleculeFromSmiles(std::declval<const std::string &>())
    } -> std::same_as<std::unique_ptr<typename FrameworkT::QueryMoleculeT>>;

    {
    framework.fingerprintFromMolecule(std::declval<const typename FrameworkT::MoleculeT &>())
    } -> std::same_as<std::unique_ptr<typename FrameworkT::FingerprintT>>;

    {
    framework.queryFingerprintFromFingerprint(std::declval<const typename FrameworkT::FingerprintT &>())
    } -> std::same_as<std::unique_ptr<typename FrameworkT::QueryFingerprintT>>;

    {
    framework.getEmptyFingerprint()
    } -> std::same_as<typename FrameworkT::FingerprintT>;

    {
    framework.compressMolecule(std::declval<const typename FrameworkT::MoleculeT &>())
    } -> std::same_as<std::unique_ptr<typename FrameworkT::StorageMoleculeT>>;

    {
    framework.decompressMolecule(std::declval<const typename FrameworkT::StorageMoleculeT &>())
    } -> std::same_as<std::unique_ptr<typename FrameworkT::MoleculeT>>;

    {
    framework.isSubstructure(std::declval<const typename FrameworkT::QueryMoleculeT &>(),
                             std::declval<const typename FrameworkT::MoleculeT &>())
    } -> std::convertible_to<bool>;

    {
    framework.isSubFingerprint(std::declval<const typename FrameworkT::QueryFingerprintT &>(),
                               std::declval<const typename FrameworkT::FingerprintT &>())
    } -> std::convertible_to<bool>;

    {
    framework.getFingerprintBit(std::declval<const typename FrameworkT::FingerprintT &>(),
                                std::declval<size_t>())
    } -> std::convertible_to<bool>;

    {
    framework.setFingerprintBit(std::declval<typename FrameworkT::FingerprintT &>(),
                                std::declval<size_t>(),
                                std::declval<bool>())
    } -> std::same_as<void>;

    {
    framework.getFingerprintSize()
    } -> std::same_as<size_t>;
};