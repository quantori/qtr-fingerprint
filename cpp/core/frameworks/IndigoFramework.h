#pragma once

#include <string>

#include "BingoNoSQL.h"
#include "IndigoSession.h"
#include "indigo.h"
#include "base_cpp/array.h"
#include "frameworks/FrameworkInterface.h"
#include "molecule/molecule.h"
#include "molecule/query_molecule.h"
#include "IndigoQueryFingerprint.h"
#include "molecule/molecule_fingerprint.h"

class IndigoFramework {
public:
    using MoleculeT = indigo::Molecule;
    using QueryMoleculeT = indigo::QueryMolecule;
    using StorageMoleculeT = std::string;

    using FingerprintInnerT = byte;
    using FingerprintT = indigo::Array<FingerprintInnerT>;
    using QueryFingerprintT = IndigoQueryFingerprint<FingerprintInnerT>;

    void init(const Config &config);

    static IndigoFramework &getInstance();

    static std::unique_ptr<MoleculeT> moleculeFromSmiles(const std::string &smiles);

    static std::string moleculeToSmiles(const MoleculeT &molecule);

    static std::unique_ptr<QueryMoleculeT> queryMoleculeFromSmiles(const std::string &smiles);

    [[nodiscard]] std::unique_ptr<FingerprintT> fingerprintFromMolecule(const MoleculeT &molecule) const;

    static std::unique_ptr<QueryFingerprintT> queryFingerprintFromFingerprint(const FingerprintT &fingerprint);

    static std::unique_ptr<StorageMoleculeT> compressMolecule(const MoleculeT &molecule);

    static std::unique_ptr<MoleculeT> decompressMolecule(const StorageMoleculeT &compressedMolecule);

    static bool isSubstructure(const QueryMoleculeT &queryMolecule, const MoleculeT &molecule);

    static bool getFingerprintBit(const FingerprintT &fingerprint, size_t idx);

    [[nodiscard]] size_t getFingerprintSize() const;

    static void setFingerprintBit(FingerprintT &fingerprint, size_t idx, bool val);

    static bool isSubFingerprint(const QueryFingerprintT &fingerprint1, const FingerprintT &fingerprint2);

    [[nodiscard]] FingerprintT getEmptyFingerprint() const;

private:
    double _fingerprintRatio = 1.0;

    [[nodiscard]] indigo::MoleculeFingerprintParameters getFingerprintParams() const;
    [[nodiscard]] std::unique_ptr<FingerprintT> tryBuildFingerprintFromMolecule(const MoleculeT &molecule) const;
};

static_assert(FrameworkInterface<IndigoFramework>, "IndigoFramework must satisfy FrameworkInterface");
