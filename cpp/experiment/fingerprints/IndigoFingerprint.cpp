#include "IndigoFingerprint.h"

namespace {
    const Indigo indigoInstance;
}

IndigoFingerprint::IndigoFingerprint(const MoleculeType &mol) : IndigoFingerprint(mol.smiles()) {
    // TODO: is it possible to avoid smiles cast?
}

IndigoFingerprint::IndigoFingerprint() {
    _fingerprint = std::make_unique<indigo::Array<byte>>();
    _fingerprint->resize(size() / 8);
}

bool IndigoFingerprint::isSubFingerprintOf(const IndigoFingerprint &other) const {
    assert(_fingerprint->size() * 8 == size());
    assert(other.size() == size());
    bool res = true;
    for (int i = 0; i < (int) _fingerprint->size() && res; i++) {
        res &= _fingerprint->operator[](i) <= other._fingerprint->operator[](i);
    }
    return res;
}

bool IndigoFingerprint::getBit(size_t index) const {
    byte elem = _fingerprint->operator[]((int) index / 8);
    return bool(elem >> (index % 8) & 1);
}

size_t IndigoFingerprint::size() {
    return 3736;
}

IndigoFingerprint &IndigoFingerprint::operator|=(const IndigoFingerprint &other) {
    assert(_fingerprint->size() == other._fingerprint->size());
    assert(_fingerprint->size() * 8 == size());
    for (int i = 0; i < _fingerprint->size(); i++) {
        _fingerprint->operator[](i) |= other._fingerprint->operator[](i);
    }
    return *this;
}

IndigoFingerprint::IndigoFingerprint(const std::string &smiles) {
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
    assert(_fingerprint->size() * 8 == size());
}
