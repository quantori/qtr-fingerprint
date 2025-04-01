#include "RDKitFingerprint.h"
#include "DataStructs/BitOps.h"

#include "GraphMol/Fingerprints/Fingerprints.h"

#include "Profiling.h"

#define RDKIT_FP_LEN (1048576)

bool RDKitFingerprint::isSubFingerprintOf(const RDKitFingerprint &other) const {
    ProfileScope("RDKitFingerprint::isSubFingerprintOf");
    return AllProbeBitsMatch(*_fingerprint, *other._fingerprint);
}

bool RDKitFingerprint::getBit(size_t index) const {
    return _fingerprint->getBit(index);
}

// RDKitFingerprint::RDKitFingerprint(const RDKit::ROMol &mol) : _fingerprint(RDKit::LayeredFingerprintMol(mol, 0xFFFFFFFF, 1, 10, RDKIT_FP_LEN)) {
RDKitFingerprint::RDKitFingerprint(const RDKit::ROMol &mol) : _fingerprint(RDKit::PatternFingerprintMol(mol, RDKIT_FP_LEN)) {
}

const ExplicitBitVect &RDKitFingerprint::bitVector() const {
    return *_fingerprint;
}

size_t RDKitFingerprint::size() const {
    return _fingerprint->size();
}

RDKitFingerprint& RDKitFingerprint::operator|=(const RDKitFingerprint &other) {
    *_fingerprint |= *other._fingerprint;
    return *this;
}

RDKitFingerprint::RDKitFingerprint() {
    _fingerprint = std::make_unique<ExplicitBitVect>(RDKIT_FP_LEN);
}
