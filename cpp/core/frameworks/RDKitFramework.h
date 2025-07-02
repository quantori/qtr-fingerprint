#pragma once

#include <string>

#include "GraphMol/ROMol.h"
#include "GraphMol/Fingerprints/Fingerprints.h"
#include "GraphMol/Substruct/SubstructMatch.h"
#include "DataStructs/ExplicitBitVect.h"

#include "frameworks/FrameworkInterface.h"
#include "Bitset.h"

class RDKitQueryFingerprint {
public:
    explicit RDKitQueryFingerprint(const ExplicitBitVect &fingerprint);

    bool isSubFingerprint(const ExplicitBitVect& fingerprint) const;

private:
    std::vector<int> _bits;
};

class RDKitFramework {
public:
    using MoleculeT = RDKit::ROMol;
    using QueryMoleculeT = RDKit::ROMol;
    using StorageMoleculeT = std::string;

    using FingerprintT = ExplicitBitVect;
    using QueryFingerprintT = RDKitQueryFingerprint;

    void init(const Config& config);

    static RDKitFramework & getInstance();

    static std::unique_ptr<MoleculeT> moleculeFromSmiles(const std::string &smiles);

    static std::string moleculeToSmiles(const MoleculeT &molecule);

    static std::unique_ptr<QueryMoleculeT> queryMoleculeFromSmiles(const std::string &smiles);

    [[nodiscard]] std::unique_ptr<FingerprintT> fingerprintFromMolecule(const MoleculeT &molecule) const;

    static std::unique_ptr<QueryFingerprintT> queryFingerprintFromFingerprint(const FingerprintT &fingerprint);

    static std::unique_ptr<StorageMoleculeT> compressMolecule(const MoleculeT &molecule);

    static std::unique_ptr<MoleculeT> decompressMolecule(const StorageMoleculeT &compressedMolecule);

    static bool isSubstructure(const QueryMoleculeT &queryMolecule, const MoleculeT &molecule);

    static RDKit::SubstructMatchParameters getSubstructMatchParameters();

    static bool getFingerprintBit(const FingerprintT &fingerprint, size_t idx);

    [[nodiscard]] size_t getFingerprintSize() const;

    static void setFingerprintBit(FingerprintT &fingerprint, size_t idx, bool val);

    static bool isSubFingerprint(const QueryFingerprintT &fingerprint1, const FingerprintT &fingerprint2);

    FingerprintT getEmptyFingerprint() const;

private:
    size_t _fingerprintLength = 2048;
};

static_assert(FrameworkInterface<RDKitFramework>, "RDKitFramework must satisfy FrameworkInterface");
