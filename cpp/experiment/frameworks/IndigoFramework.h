#pragma once

#include <string>

#include "BingoNoSQL.h"
#include "IndigoSession.h"
#include "indigo.h"
#include "base_cpp/array.h"

#include "frameworks/FrameworkInterface.h"

class IndigoFramework {
public:
    using FingerprintT = indigo::Array<byte>;
    using MoleculeT = indigo_cpp::IndigoMolecule;
    using StorageMoleculeT = indigo_cpp::IndigoMolecule;
    using QueryMoleculeT = indigo_cpp::IndigoQueryMolecule;

    static std::unique_ptr<MoleculeT> moleculeFromSmiles(const std::string &smiles);

    static std::string moleculeToSmiles(const MoleculeT &molecule);

    static std::unique_ptr<QueryMoleculeT> queryMoleculeFromSmiles(const std::string &smiles);

    static std::unique_ptr<FingerprintT> fingerprintFromMolecule(const MoleculeT &molecule);

    static std::unique_ptr<StorageMoleculeT> compressMolecule(const MoleculeT &molecule);

    static std::unique_ptr<MoleculeT> decompressMolecule(const StorageMoleculeT &compressedMolecule);

    static bool isSubstructure(const QueryMoleculeT &queryMolecule, const MoleculeT &molecule);

    static bool getFingerprintBit(const FingerprintT &fingerprint, size_t idx);

    static size_t getFingerprintSize();

    static void setFingerprintBit(FingerprintT &fingerprint, size_t idx, bool val);

    static bool isSubFingerprint(const FingerprintT &fingerprint1, const FingerprintT &fingerprint2);

    static std::shared_ptr<indigo_cpp::IndigoSession> getGlobalIndigoSession();

    static FingerprintT getEmptyFingerprint();
};

static_assert(FrameworkInterface<IndigoFramework>, "IndigoFramework must satisfy FrameworkInterface");
