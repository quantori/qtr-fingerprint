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
RDKitFramework::fingerprintFromMolecule(const MoleculeT &molecule) const {
    return std::unique_ptr<ExplicitBitVect>(
            RDKit::PatternFingerprintMol(molecule, getFingerprintSize())
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

size_t RDKitFramework::getFingerprintSize() const {
    return _fingerprintLength;
}

void RDKitFramework::setFingerprintBit(RDKitFramework::FingerprintT &fingerprint, size_t idx, bool val) {
    if (val) {
        fingerprint.setBit(idx);
    } else {
        fingerprint.unsetBit(idx);
    }
}

bool RDKitFramework::isSubFingerprint(const RDKitFramework::QueryFingerprintT &fingerprint1,
                                      const RDKitFramework::FingerprintT &fingerprint2) {
    ProfileScope("RDKitFramework::isSubFingerprint");
    return fingerprint1.isSubFingerprint(fingerprint2);
}

std::string RDKitFramework::moleculeToSmiles(const RDKitFramework::MoleculeT &molecule) {
    return RDKit::MolToSmiles(molecule);
}

RDKitFramework::FingerprintT RDKitFramework::getEmptyFingerprint() const {
    return FingerprintT(getFingerprintSize());
}

std::unique_ptr<RDKitFramework::QueryFingerprintT>
RDKitFramework::queryFingerprintFromFingerprint(const RDKitFramework::FingerprintT &fingerprint) {
    return std::make_unique<QueryFingerprintT>(fingerprint);
}

void RDKitFramework::init(const Config &config) {
    _fingerprintLength = size_t(ceil(2048 * config.getDouble("fpRatio", 1)));
}

RDKitFramework & RDKitFramework::getInstance() {
    static RDKitFramework instance;
    return instance;
}

RDKitQueryFingerprint::RDKitQueryFingerprint(const ExplicitBitVect &fingerprint) {
    fingerprint.getOnBits(_bits);
}

bool RDKitQueryFingerprint::isSubFingerprint(const ExplicitBitVect &fingerprint) const {
    for (auto &i: _bits) {
        if (!fingerprint.getBit(i)) {
            return false;
        }
    }
    return true;
}
