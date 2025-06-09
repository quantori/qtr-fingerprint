#pragma once

#include <string>

#include "BingoNoSQL.h"
#include "IndigoSession.h"
#include "indigo.h"
#include "base_cpp/array.h"

#include "frameworks/FrameworkInterface.h"
#include "molecule/molecule.h"
#include "IndigoQueryFingerprint.h"

class IndigoCppFramework {
public:
    using MoleculeT = indigo_cpp::IndigoMolecule;
    using QueryMoleculeT = indigo_cpp::IndigoQueryMolecule;
    using StorageMoleculeT = indigo_cpp::IndigoMolecule;

    using FingerprintT = indigo::Array<byte>;
    using QueryFingerprintT = IndigoQueryFingerprint<byte>;

    void init(const Config& config);

    static IndigoCppFramework& getInstance();

    std::unique_ptr<MoleculeT> moleculeFromSmiles(const std::string &smiles);

    static std::string moleculeToSmiles(const MoleculeT &molecule);

    [[nodiscard]] std::unique_ptr<QueryMoleculeT> queryMoleculeFromSmiles(const std::string &smiles) const;

    std::unique_ptr<FingerprintT> fingerprintFromMolecule(const MoleculeT &molecule);

    static std::unique_ptr<QueryFingerprintT> queryFingerprintFromFingerprint(const FingerprintT &fingerprint);

    static std::unique_ptr<StorageMoleculeT> compressMolecule(const MoleculeT &molecule);

    static std::unique_ptr<MoleculeT> decompressMolecule(const StorageMoleculeT &compressedMolecule);

    static bool isSubstructure(const QueryMoleculeT &queryMolecule, const MoleculeT &molecule);

    static bool getFingerprintBit(const FingerprintT &fingerprint, size_t idx);

    static size_t getFingerprintSize();

    static void setFingerprintBit(FingerprintT &fingerprint, size_t idx, bool val);

    static bool isSubFingerprint(const QueryFingerprintT &fingerprint1, const FingerprintT &fingerprint2);

    static FingerprintT getEmptyFingerprint();

    std::shared_ptr<indigo_cpp::IndigoSession> getSession();

private:
    static std::unique_ptr<IndigoCppFramework::FingerprintT>
    tryBuildFingerprintFromMolecule(const IndigoCppFramework::MoleculeT &molecule);
};

static_assert(FrameworkInterface<IndigoCppFramework>, "IndigoCppFramework must satisfy FrameworkInterface");
