#pragma once

#include <ranges>
#include <execution>

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
        auto range = std::views::iota(size_t(0), smiles.size());
        std::mutex mutex;
        tbb::parallel_for(
            tbb::blocked_range<size_t>(range.begin(), range.end()),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t idx = r.begin(); idx != r.end(); ++idx) {
                    auto& s = smiles.smiles(idx);
                    auto molecule = FrameworkT::moleculeFromSmiles(s);
                    auto compressedMol = FrameworkT::compressMolecule(*molecule);
                    auto fingerprint = FrameworkT::fingerprintFromMolecule(*molecule);
                    {
                        std::lock_guard<std::mutex> lockGuard(mutex);
                        _fingerprints.push_back(std::move(fingerprint));
                        _storedMolecules.push_back(std::move(compressedMol));
                    }
                }
        });
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
