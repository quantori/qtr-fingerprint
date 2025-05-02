#pragma once

#include <concepts>
#include <string>
#include <memory>


template<typename FrameworkT>
concept FrameworkInterface = requires {
    typename FrameworkT::MoleculeT;
    typename FrameworkT::QueryMoleculeT;
    typename FrameworkT::StorageMoleculeT;

    typename FrameworkT::FingerprintT;
    typename FrameworkT::QueryFingerprintT;

//    requires std::is_copy_constructible_v<typename FrameworkT::MoleculeT>;
//    requires std::is_copy_constructible_v<typename FrameworkT::QueryMoleculeT>;
//    requires std::is_copy_constructible_v<typename FrameworkT::StorageMoleculeT>;

    {
    FrameworkT::moleculeFromSmiles(std::declval<const std::string &>())
    } -> std::same_as<std::unique_ptr<typename FrameworkT::MoleculeT>>;

    {
    FrameworkT::moleculeToSmiles(std::declval<typename FrameworkT::MoleculeT &>())
    } -> std::convertible_to<std::string>;

    {
    FrameworkT::queryMoleculeFromSmiles(std::declval<const std::string &>())
    } -> std::same_as<std::unique_ptr<typename FrameworkT::QueryMoleculeT>>;

    {
    FrameworkT::fingerprintFromMolecule(std::declval<const typename FrameworkT::MoleculeT &>())
    } -> std::same_as<std::unique_ptr<typename FrameworkT::FingerprintT>>;

    {
    FrameworkT::queryFingerprintFromFingerprint(std::declval<const typename FrameworkT::FingerprintT &>())
    } -> std::same_as<std::unique_ptr<typename FrameworkT::QueryFingerprintT>>;

    {
    FrameworkT::getEmptyFingerprint()
    } -> std::same_as<typename FrameworkT::FingerprintT>;

    {
    FrameworkT::compressMolecule(std::declval<const typename FrameworkT::MoleculeT &>())
    } -> std::same_as<std::unique_ptr<typename FrameworkT::StorageMoleculeT>>;

    {
    FrameworkT::decompressMolecule(std::declval<const typename FrameworkT::StorageMoleculeT &>())
    } -> std::same_as<std::unique_ptr<typename FrameworkT::MoleculeT>>;

    {
    FrameworkT::isSubstructure(std::declval<const typename FrameworkT::QueryMoleculeT &>(),
                               std::declval<const typename FrameworkT::MoleculeT &>())
    } -> std::convertible_to<bool>;

    {
    FrameworkT::isSubFingerprint(std::declval<const typename FrameworkT::QueryFingerprintT &>(),
                                 std::declval<const typename FrameworkT::FingerprintT &>())
    } -> std::convertible_to<bool>;

    {
    FrameworkT::getFingerprintBit(std::declval<const typename FrameworkT::FingerprintT &>(),
                                  std::declval<size_t>())
    } -> std::convertible_to<bool>;

    {
    FrameworkT::setFingerprintBit(std::declval<typename FrameworkT::FingerprintT &>(),
                                  std::declval<size_t>(),
                                  std::declval<bool>())
    } -> std::same_as<void>;

    {
    FrameworkT::getFingerprintSize()
    } -> std::same_as<size_t>;
};