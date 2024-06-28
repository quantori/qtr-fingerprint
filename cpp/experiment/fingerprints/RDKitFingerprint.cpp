#include "RDKitFingerprint.h"
#include "DataStructs/BitOps.h"

#include "GraphMol/Fingerprints/Fingerprints.h"

bool RDKitFingerprint::isSubFingerprintOf(const RDKitFingerprint &other) const {
    return AllProbeBitsMatch(*_fingerprint, *other._fingerprint);
}

bool RDKitFingerprint::getBit(size_t index) const {
    return _fingerprint->getBit(index);
}

RDKitFingerprint::RDKitFingerprint(const RDKit::ROMol &mol) : _fingerprint(RDKit::PatternFingerprintMol(mol)) {
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
    _fingerprint = std::make_unique<ExplicitBitVect>(2048);
}
