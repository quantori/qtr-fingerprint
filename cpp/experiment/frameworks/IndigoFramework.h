#pragma once

#include <string>

#include "BingoNoSQL.h"
#include "IndigoSession.h"
#include "indigo.h"
#include "base_cpp/array.h"

#include "frameworks/FrameworkInterface.h"
#include "molecule/molecule.h"
#include "IndigoQueryFingerprint.h"

class IndigoFramework {
public:
    using MoleculeT = indigo::Molecule;
    using QueryMoleculeT = indigo::QueryMolecule;
    using StorageMoleculeT = std::string;

    using FingerprintInnerT = byte;
    using FingerprintT = indigo::Array<FingerprintInnerT>;
    using QueryFingerprintT = IndigoQueryFingerprint<FingerprintInnerT>;

    explicit IndigoFramework(const Config& config);

    static std::unique_ptr<MoleculeT> moleculeFromSmiles(const std::string &smiles);

    static std::string moleculeToSmiles(const MoleculeT &molecule);

    static std::unique_ptr<QueryMoleculeT> queryMoleculeFromSmiles(const std::string &smiles);

    static std::unique_ptr<FingerprintT> fingerprintFromMolecule(const MoleculeT &molecule);

    static std::unique_ptr<QueryFingerprintT> queryFingerprintFromFingerprint(const FingerprintT& fingerprint);

    static std::unique_ptr<StorageMoleculeT> compressMolecule(const MoleculeT &molecule);

    static std::unique_ptr<MoleculeT> decompressMolecule(const StorageMoleculeT &compressedMolecule);

    static bool isSubstructure(const QueryMoleculeT &queryMolecule, const MoleculeT &molecule);

    static bool getFingerprintBit(const FingerprintT &fingerprint, size_t idx);

    static size_t getFingerprintSize();

    static void setFingerprintBit(FingerprintT &fingerprint, size_t idx, bool val);

    static bool isSubFingerprint(const QueryFingerprintT &fingerprint1, const FingerprintT &fingerprint2);

    static FingerprintT getEmptyFingerprint();
};

static_assert(FrameworkInterface<IndigoFramework>, "IndigoFramework must satisfy FrameworkInterface");
