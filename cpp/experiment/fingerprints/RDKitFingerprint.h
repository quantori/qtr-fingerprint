#pragma once

#include "GraphMol/Fingerprints/Fingerprints.h"
#include "FingerprintConcept.h"

class RDKitFingerprint {
public:
    using MoleculeType = RDKit::ROMol;

    explicit RDKitFingerprint(const MoleculeType &mol);

    RDKitFingerprint();

    [[nodiscard]] bool isSubFingerprintOf(const RDKitFingerprint &other) const;

    [[nodiscard]] bool getBit(size_t index) const;

    [[nodiscard]] const ExplicitBitVect &bitVector() const;

    [[nodiscard]] size_t size() const;

    RDKitFingerprint &operator|=(const RDKitFingerprint &other);

private:
    std::unique_ptr<ExplicitBitVect> _fingerprint;
};
