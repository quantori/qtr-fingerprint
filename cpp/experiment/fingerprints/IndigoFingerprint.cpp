#include "IndigoFingerprint.h"

namespace {
    const Indigo indigoInstance;
}

IndigoFingerprint::IndigoFingerprint(indigo_cpp::IndigoMolecule &mol) {
    // TODO: is it possible to cast indigo_cpp::IndigoMolecule to indigo::Molecule directly?
    std::string smiles = mol.smiles();
    indigo::BufferScanner scanner(smiles.c_str(), smiles.size(), false);
    indigo::SmilesLoader loader(scanner);
    indigo::Molecule molecule;
    molecule.aromatize(indigo::AromaticityOptions());
    assert(molecule.isAromatized());
    bingo::IndexMolecule indexMolecule(molecule, indigo::AromaticityOptions());
    _fingerprint = std::make_unique<indigo::Array<byte>>();
    indigo::MoleculeFingerprintBuilder fingerprintBuilder(molecule, indigoInstance.fp_params);
    fingerprintBuilder.parseFingerprintType("sub", false);
    fingerprintBuilder.process();
    _fingerprint->copy(fingerprintBuilder.get(), indigoInstance.fp_params.fingerprintSize());
    assert(_fingerprint->size() == size());
}

IndigoFingerprint::IndigoFingerprint() {
    _fingerprint = std::make_unique<indigo::Array<byte>>();
    _fingerprint->resize(size());
}

bool IndigoFingerprint::isSubFingerprintOf(const IndigoFingerprint &other) const {
    assert(_fingerprint->size() == size());
    assert(other.size() == size());
    bool res = true;
    for (int i = 0; i < (int)IndigoFingerprint::size() && res; i++) {
        res &= _fingerprint->operator[](i) <=  other._fingerprint->operator[](i);
    }
    return res;
}

bool IndigoFingerprint::getBit(size_t index) const {
    byte elem = _fingerprint->operator[]((int)index / 8);
    return bool(elem >> (index % 8) & 1);
}

size_t IndigoFingerprint::size() {
    return 3736;
}

IndigoFingerprint &IndigoFingerprint::operator|=(const IndigoFingerprint &other) {
    for (int i = 0; i < (int)IndigoFingerprint::size(); i++) {
        _fingerprint->operator[](i) |= other._fingerprint->operator[](i);
    }
    return *this;
}
