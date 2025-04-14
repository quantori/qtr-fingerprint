#include "RDKitFramework.h"

#include <glog/logging.h>

#include "GraphMol/SubstructLibrary/SubstructLibrary.h"
#include "Profiling.h"

std::unique_ptr<RDKitFramework::MoleculeT> RDKitFramework::moleculeFromSmiles(const std::string &smiles) {
    try {
        auto res = std::unique_ptr<RDKitFramework::MoleculeT>(RDKit::SmilesToMol(smiles));
        return res;
    } catch (const std::exception &e) {
        LOG(WARNING) << "Cannot parse SMILES: " << smiles << " error: " << e.what();
        return nullptr;
    }
}

std::unique_ptr<RDKitFramework::QueryMoleculeT> RDKitFramework::queryMoleculeFromSmiles(const std::string &smiles) {
    return moleculeFromSmiles(smiles);
}

std::unique_ptr<RDKitFramework::FingerprintT>
RDKitFramework::fingerprintFromMolecule(const RDKitFramework::MoleculeT &molecule) {
    return std::unique_ptr<RDKitFramework::FingerprintT>(
            RDKit::PatternFingerprintMol(molecule, RDKitFramework::getFingerprintSize())
    );
}

std::unique_ptr<RDKitFramework::StorageMoleculeT>
RDKitFramework::compressMolecule(const MoleculeT &molecule) {
    auto pickle = std::make_unique<std::string>();
    RDKit::MolPickler::pickleMol(molecule, *pickle);
    return pickle;
}

std::unique_ptr<RDKitFramework::MoleculeT>
RDKitFramework::decompressMolecule(const RDKitFramework::StorageMoleculeT &compressedMolecule) {
    ProfileScope("RDKitFramework::decompressMolecule");
    std::unique_ptr<RDKitFramework::MoleculeT> mol(new RDKit::ROMol);
    RDKit::MolPickler::molFromPickle(compressedMolecule, mol.get());
    return mol;
}

bool RDKitFramework::isSubstructure(const RDKitFramework::QueryMoleculeT &queryMolecule,
                                    const RDKitFramework::MoleculeT &molecule) {
    ProfileScope("RDKitFramework::isSubstructure");
    auto params = RDKitFramework::getSubstructMatchParameters();
    return !SubstructMatch(molecule, queryMolecule, params).empty();
}

RDKit::SubstructMatchParameters RDKitFramework::getSubstructMatchParameters() {
    RDKit::SubstructMatchParameters params;
    params.recursionPossible = true;
    params.useChirality = true;
    params.useQueryQueryMatches = false;
    return params;
}

bool RDKitFramework::getFingerprintBit(const RDKitFramework::FingerprintT &fingerprint, size_t idx) {
    return fingerprint.getBit(idx);
}

size_t RDKitFramework::getFingerprintSize() {
    return 2048;
}

void RDKitFramework::setFingerprintBit(RDKitFramework::FingerprintT &fingerprint, size_t idx, bool val) {
    if (val) {
        fingerprint.setBit(idx);
    } else {
        fingerprint.unsetBit(idx);
    }
}

bool RDKitFramework::isSubFingerprint(const RDKitFramework::FingerprintT &fingerprint1,
                                      const RDKitFramework::FingerprintT &fingerprint2) {
    ProfileScope("RDKitFramework::isSubFingerprint");
    for (size_t i = 0; i < RDKitFramework::getFingerprintSize(); i++) {
        if (RDKitFramework::getFingerprintBit(fingerprint1, i) > RDKitFramework::getFingerprintBit(fingerprint2, i)) {
            return false;
        }
    }
    return true;
}

std::string RDKitFramework::moleculeToSmiles(const RDKitFramework::MoleculeT &molecule) {
    return RDKit::MolToSmiles(molecule);
}
