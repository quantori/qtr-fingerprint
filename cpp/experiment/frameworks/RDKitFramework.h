#pragma once

#include <string>

#include "GraphMol/ROMol.h"
#include "GraphMol/Fingerprints/Fingerprints.h"
#include "GraphMol/Substruct/SubstructMatch.h"

#include "frameworks/FrameworkInterface.h"

class RDKitFramework {
public:
    using MoleculeT = RDKit::ROMol;
    using FingerprintT = ExplicitBitVect;
    using QueryMoleculeT = RDKit::ROMol;
    using StorageMoleculeT = std::string;

    static std::unique_ptr<MoleculeT> moleculeFromSmiles(const std::string &smiles);

    static std::string moleculeToSmiles(const MoleculeT &molecule);

    static std::unique_ptr<QueryMoleculeT> queryMoleculeFromSmiles(const std::string &smiles);

    static std::unique_ptr<FingerprintT> fingerprintFromMolecule(const MoleculeT &molecule);

    static std::unique_ptr<StorageMoleculeT> compressMolecule(const MoleculeT &molecule);

    static std::unique_ptr<MoleculeT> decompressMolecule(const StorageMoleculeT &compressedMolecule);

    static bool isSubstructure(const QueryMoleculeT &queryMolecule, const MoleculeT &molecule);

    static RDKit::SubstructMatchParameters getSubstructMatchParameters();

    static bool getFingerprintBit(const FingerprintT& fingerprint, size_t idx);

    static size_t getFingerprintSize();

    static void setFingerprintBit(FingerprintT& fingerprint, size_t idx, bool val);

    static bool isSubFingerprint(const FingerprintT& fingerprint1, const FingerprintT& fingerprint2);
};

static_assert(FrameworkInterface<RDKitFramework>, "RDKitFramework must satisfy FrameworkInterface");
