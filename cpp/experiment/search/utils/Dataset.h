#pragma once

#include <vector>
#include <memory>
#include <execution>
#include <ranges>

#include <glog/logging.h>

#include "frameworks/FrameworkInterface.h"

template<typename FrameworkT> requires FrameworkInterface<FrameworkT>
class Dataset {
public:
    using MoleculeT = typename FrameworkT::MoleculeT;
    using StorageMoleculeT = typename FrameworkT::StorageMoleculeT;
    using FingerprintT = typename FrameworkT::FingerprintT;

    explicit Dataset(const std::vector<std::string> &smilesSequence) : _molecules() {
        std::mutex mutex;
        std::for_each(std::execution::par, smilesSequence.begin(), smilesSequence.end(), [&](const std::string &s) {
            auto mol = FrameworkT::moleculeFromSmiles(s);
            if (mol == nullptr) {
                LOG(WARNING) << "Skip unparsed SMILES: " << s;
                return;
            }
            {
                std::lock_guard<std::mutex> lockGuard(mutex);
                _molecules.push_back(mol);
            }
        });
    }

    void buildFingerprints() {
        _fingerprints.resize(size());
        std::for_each(std::execution::par, std::views::iota(0, size()), [&](int idx) {
            auto &mol = *_molecules[idx];
            _fingerprints[idx] = FrameworkT::fingerprintFromMolecule(mol);
        });
    }

    [[nodiscard]] size_t size() const {
        return _molecules.size();
    }

    const FingerprintT& fingerprint(size_t idx) {
        return *_fingerprints.at(idx);
    }

    const StorageMoleculeT& storageMolecule(size_t idx) {
        return *_molecules.at(idx);
    }

    std::unique_ptr<MoleculeT> molecule(size_t idx) {
        const auto& storageMol = storageMolecule(idx);
        auto mol = FrameworkT::decompressMolecule(storageMol);
        return mol;
    }

private:
    std::vector<std::unique_ptr<StorageMoleculeT>> _molecules;
    std::vector<std::unique_ptr<FingerprintT>> _fingerprints;
};