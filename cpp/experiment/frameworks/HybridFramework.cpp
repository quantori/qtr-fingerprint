#include "frameworks/HybridFramework.h"

#include "IndigoSession.h"
#include "indigo.h"
#include "base_cpp/array.h"
#include "molecule/query_molecule.h"

#include "GraphMol/ROMol.h"
#include "GraphMol/Fingerprints/Fingerprints.h"
#include "GraphMol/SmilesParse/SmilesParse.h"
#include "GraphMol/SmilesParse/SmilesWrite.h"
#include "GraphMol/Substruct/SubstructMatch.h"

// Molecule operations use Indigo
std::unique_ptr<HybridFramework::MoleculeT> HybridFramework::moleculeFromSmiles(const std::string &smiles) {
    return IndigoFramework::moleculeFromSmiles(smiles);
}

std::string HybridFramework::moleculeToSmiles(const MoleculeT &molecule) {
    return IndigoFramework::moleculeToSmiles(molecule);
}

std::unique_ptr<HybridFramework::QueryMoleculeT> HybridFramework::queryMoleculeFromSmiles(const std::string &smiles) {
    return IndigoFramework::queryMoleculeFromSmiles(smiles);
}

std::unique_ptr<HybridFramework::StorageMoleculeT> HybridFramework::compressMolecule(const MoleculeT &molecule) {
    return IndigoFramework::compressMolecule(molecule);
}

std::unique_ptr<HybridFramework::MoleculeT> HybridFramework::decompressMolecule(const StorageMoleculeT &compressedMolecule) {
    return IndigoFramework::decompressMolecule(compressedMolecule);
}

bool HybridFramework::isSubstructure(const QueryMoleculeT &queryMolecule, const MoleculeT &molecule) {
    return IndigoFramework::isSubstructure(queryMolecule, molecule);
}

// Fingerprint operations use RDKit
std::unique_ptr<HybridFramework::FingerprintT> HybridFramework::fingerprintFromMolecule(const MoleculeT &molecule) {
    // Convert Indigo molecule to SMILES
    std::string smiles = moleculeToSmiles(molecule);
    
    // Convert SMILES to RDKit molecule
    std::unique_ptr<RDKit::ROMol> rdkitMol = RDKitFramework::moleculeFromSmiles(smiles);
    
    // Generate fingerprint using RDKit
    return RDKitFramework::fingerprintFromMolecule(*rdkitMol);
}

std::unique_ptr<HybridFramework::QueryFingerprintT> HybridFramework::queryFingerprintFromFingerprint(const FingerprintT &fingerprint) {
    return RDKitFramework::queryFingerprintFromFingerprint(fingerprint);
}

bool HybridFramework::getFingerprintBit(const FingerprintT &fingerprint, size_t idx) {
    return RDKitFramework::getFingerprintBit(fingerprint, idx);
}

size_t HybridFramework::getFingerprintSize() {
    return RDKitFramework::getFingerprintSize();
}

void HybridFramework::setFingerprintBit(FingerprintT &fingerprint, size_t idx, bool val) {
    RDKitFramework::setFingerprintBit(fingerprint, idx, val);
}

bool HybridFramework::isSubFingerprint(const QueryFingerprintT &fingerprint1, const FingerprintT &fingerprint2) {
    return RDKitFramework::isSubFingerprint(fingerprint1, fingerprint2);
}

HybridFramework::FingerprintT HybridFramework::getEmptyFingerprint() {
    return RDKitFramework::getEmptyFingerprint();
} 