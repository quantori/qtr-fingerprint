#include "IndigoFingerprint.h"

#include "Profiling.h"

namespace {
    const Indigo indigoInstance;
}

IndigoFingerprint::IndigoFingerprint(const MoleculeType &mol) {
    // Maybe one needs to aromatize
    //    mol.aromatize(indigo::AromaticityOptions());
    indigo::MoleculeFingerprintBuilder fingerprintBuilder((MoleculeType &)mol, indigoInstance.fp_params);
    fingerprintBuilder.parseFingerprintType("sub", false);
    fingerprintBuilder.process();
    indigo::Array<byte> rawFingerprint;
    rawFingerprint.copy(fingerprintBuilder.get(), indigoInstance.fp_params.fingerprintSize());
    
    // Pack the raw fingerprint bytes into the vector of type T
    _fingerprint.resize((size() + sizeof(T) * 8 - 1) / (sizeof(T) * 8));
    for (size_t i = 0; i < rawFingerprint.size(); ++i) {
        _fingerprint[i / sizeof(T)] |= static_cast<T>(rawFingerprint[i]) << ((i % sizeof(T)) * 8);
    }
}

IndigoFingerprint::IndigoFingerprint() {
    _fingerprint.resize((size() + sizeof(T) * 8 - 1) / (sizeof(T) * 8));
}

bool IndigoFingerprint::isSubFingerprintOf(const IndigoFingerprint &other) const {
    ProfileScope("IndigoFingerprint::isSubFingerprintOf");
    assert(_fingerprint.size() == other._fingerprint.size());
    
    for (size_t i = 0; i < _fingerprint.size(); ++i) {
        if ((_fingerprint[i] & other._fingerprint[i]) != _fingerprint[i]) {
            return false;
        }
    }
    return true;
}

bool IndigoFingerprint::getBit(size_t index) const {
    return bool(_fingerprint[index / (sizeof(T) * 8)] >> (index % (sizeof(T) * 8)) & 1);
}

size_t IndigoFingerprint::size() {
    return 3736;
}

IndigoFingerprint &IndigoFingerprint::operator|=(const IndigoFingerprint &other) {
    assert(_fingerprint.size() == other._fingerprint.size());
    assert(_fingerprint.size() * 8 == size());
    for (int i = 0; i < _fingerprint.size(); i++) {
        _fingerprint[i] |= other._fingerprint[i];
    }
    return *this;
}

IndigoFingerprint::IndigoFingerprint(const std::string &smiles) {
    indigo::BufferScanner scanner(smiles.c_str(), smiles.size(), false);
    indigo::SmilesLoader loader(scanner);
    indigo::Molecule molecule;
    molecule.aromatize(indigo::AromaticityOptions());
    assert(molecule.isAromatized());
    indigo::MoleculeFingerprintBuilder fingerprintBuilder(molecule, indigoInstance.fp_params);
    fingerprintBuilder.parseFingerprintType("sub", false);
    fingerprintBuilder.process();
    indigo::Array<byte> rawFingerprint;
    rawFingerprint.copy(fingerprintBuilder.get(), indigoInstance.fp_params.fingerprintSize());
    
    // Pack the raw fingerprint bytes into the vector of type T
    _fingerprint.resize((size() + sizeof(T) * 8 - 1) / (sizeof(T) * 8));
    for (size_t i = 0; i < rawFingerprint.size(); ++i) {
        _fingerprint[i / sizeof(T)] |= static_cast<T>(rawFingerprint[i]) << ((i % sizeof(T)) * 8);
    }
}
