#include "IndigoFramework.h"
#include "molecule/molecule_fingerprint.h"
#include "indigo_internal.h"
#include "IndigoSubstructureMatcher.h"
#include "molecule/molecule_substructure_matcher.h"
#include "base_cpp/scanner.h"
#include "Profiling.h"

namespace {
    const std::shared_ptr<indigo_cpp::IndigoSession> globalIndigoSession = indigo_cpp::IndigoSession::create();

    std::unique_ptr<IndigoFramework::FingerprintT>
    tryBuildFingerprintFromMolecule(const IndigoFramework::MoleculeT &molecule) {
        INDIGO_BEGIN
                {
                    indigo::Molecule &mol = self.getObject(molecule.id()).getMolecule();
                    mol.aromatize(indigo::AromaticityOptions());
                    assert(mol.isAromatized());
                    indigo::MoleculeFingerprintBuilder fingerprintBuilder(mol, self.fp_params);
                    fingerprintBuilder.parseFingerprintType("sub", false);
                    fingerprintBuilder.process();
                    auto result = std::make_unique<IndigoFramework::FingerprintT>();
                    result->reserve(self.fp_params.fingerprintSize());
                    result->copy(fingerprintBuilder.get(), self.fp_params.fingerprintSize());
                    assert(result->size() * CHAR_BIT == IndigoFramework::getFingerprintSize());
                    return result;
                }
        INDIGO_END(nullptr);
    }

}

std::unique_ptr<IndigoFramework::MoleculeT> IndigoFramework::moleculeFromSmiles(const std::string &smiles) {

    auto mol = std::make_unique<MoleculeT>(globalIndigoSession->loadMolecule(smiles));
    mol->aromatize();
    return mol;
}

std::string IndigoFramework::moleculeToSmiles(const IndigoFramework::MoleculeT &molecule) {
    return molecule.smiles();
}

std::unique_ptr<IndigoFramework::QueryMoleculeT> IndigoFramework::queryMoleculeFromSmiles(const std::string &smiles) {
    auto mol = std::make_unique<QueryMoleculeT>(globalIndigoSession->loadQueryMolecule(smiles));
    mol->aromatize();
    return mol;
}

std::unique_ptr<IndigoFramework::FingerprintT>
IndigoFramework::fingerprintFromMolecule(const IndigoFramework::MoleculeT &molecule) {
    auto res = tryBuildFingerprintFromMolecule(molecule);
    assert(res != nullptr);
    return res;
}

std::unique_ptr<IndigoFramework::StorageMoleculeT>
IndigoFramework::compressMolecule(const IndigoFramework::MoleculeT &molecule) {
    return std::make_unique<StorageMoleculeT>(molecule);
}

std::unique_ptr<IndigoFramework::MoleculeT>
IndigoFramework::decompressMolecule(const IndigoFramework::StorageMoleculeT &compressedMolecule) {
    return std::make_unique<MoleculeT>(compressedMolecule);
}

bool IndigoFramework::isSubstructure(const IndigoFramework::QueryMoleculeT &queryMolecule,
                                     const IndigoFramework::MoleculeT &molecule) {
    ProfileScope("IndigoFramework::isSubstructure");
    auto matcher = IndigoFramework::getGlobalIndigoSession()->substructureMatcher(molecule);
    return matcher.match(queryMolecule);
}

size_t IndigoFramework::getFingerprintSize() {
    return 3736;
}

bool IndigoFramework::getFingerprintBit(const IndigoFramework::FingerprintT &fingerprint, size_t idx) {
    int byteIdx = int(idx / CHAR_BIT);
    int bitIndex = int(idx % CHAR_BIT);
    unsigned char byte = fingerprint.at(byteIdx);
    return bool((byte >> bitIndex) & 1u);
}

void IndigoFramework::setFingerprintBit(IndigoFramework::FingerprintT &fingerprint, size_t idx, bool val) {
    int byteIdx = int(idx / CHAR_BIT);
    int bitIndex = int(idx % CHAR_BIT);
    unsigned char &byte = fingerprint[byteIdx];

    if (val)
        byte |= (1 << bitIndex);
    else
        byte &= ~(1 << bitIndex);
}

bool IndigoFramework::isSubFingerprint(const IndigoFramework::FingerprintT &fingerprint1,
                                       const IndigoFramework::FingerprintT &fingerprint2) {
    assert(fingerprint1.size() == fingerprint2.size());
    ProfileScope("IndigoFramework::isSubFingerprint");
    for (int i = 0; i < fingerprint1.size(); i++) {
        unsigned char c1 = fingerprint1.at(i);
        unsigned char c2 = fingerprint2.at(i);
        if ((c1 & c2) != c1) {
            return false;
        }
    }
    return true;
}

std::shared_ptr<indigo_cpp::IndigoSession> IndigoFramework::getGlobalIndigoSession() {
    return globalIndigoSession;
}

IndigoFramework::FingerprintT IndigoFramework::getEmptyFingerprint() {
    FingerprintT fp;
    int byteSize = (int) getFingerprintSize() / CHAR_BIT;
    fp.reserve(byteSize);
    fp.resize(byteSize);
    fp.zerofill();
    return fp;
}


