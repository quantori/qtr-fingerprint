#pragma once

#include "indigo.h"
#include "IndigoMolecule.h"
#include "IndigoSession.h"
#include "IndigoWriteBuffer.h"

#include "base_cpp/array.h"
#include "base_cpp/scanner.h"
#include "src/bingo_object.h"
#include "src/indigo_internal.h"
#include "molecule/smiles_loader.h"
#include "molecule/molecule_fingerprint.h"

class IndigoFingerprint {
public:
    using MoleculeType = indigo_cpp::IndigoBaseMolecule;

    explicit IndigoFingerprint(const MoleculeType &mol);

    explicit IndigoFingerprint(const std::string &smiles);

    IndigoFingerprint();

    [[nodiscard]] bool isSubFingerprintOf(const IndigoFingerprint &other) const;

    [[nodiscard]] bool getBit(size_t index) const;

    [[nodiscard]] static size_t size();

    IndigoFingerprint &operator|=(const IndigoFingerprint &other);

private:
    std::unique_ptr<indigo::Array<byte>> _fingerprint;
};
