#pragma once

#include <ranges>
#include <execution>
#include <tbb/parallel_for.h>

#include "dataset/DatasetInterface.h"
#include "frameworks/FrameworkInterface.h"

template<typename FrameworkType> requires FrameworkInterface<FrameworkType>
class CachedDataset {
public:
    using FrameworkT = FrameworkType;
    using MoleculeT = typename FrameworkT::MoleculeT;
    using StorageMoleculeT = typename FrameworkT::StorageMoleculeT;
    using FingerprintT = typename FrameworkT::FingerprintT;

    explicit CachedDataset(SmilesStorage &&smiles) {
        std::mutex mutex;
        auto &framework = FrameworkT::getInstance();
        for (size_t idx = 0; idx != smiles.size(); ++idx) {
            auto &s = smiles.smiles(idx);
            try {
                auto molecule = FrameworkT::moleculeFromSmiles(s);
                if (molecule == nullptr) {
                    continue;
                }
                auto compressedMol = framework.compressMolecule(*molecule);
                auto fingerprint = framework.fingerprintFromMolecule(*molecule);
                {
                    std::lock_guard<std::mutex> lockGuard(mutex);
                    _fingerprints.push_back(std::move(fingerprint));
                    _storedMolecules.push_back(std::move(compressedMol));
                }
            } catch (const std::exception &e) {
                LOG(ERROR) << "Error processing smiles " << s << ": " << e.what();
            }
        }
    }

    [[nodiscard]] size_t size() const {
        return _storedMolecules.size();
    }

    std::unique_ptr<MoleculeT> molecule(size_t idx) const {
        return FrameworkT::decompressMolecule(*_storedMolecules.at(idx));
    }

    const FingerprintT &fingerprint(size_t idx) const {
        return *_fingerprints.at(idx);
    }

private:
    std::vector<std::unique_ptr<StorageMoleculeT>> _storedMolecules;
    std::vector<std::unique_ptr<FingerprintT>> _fingerprints;
};
