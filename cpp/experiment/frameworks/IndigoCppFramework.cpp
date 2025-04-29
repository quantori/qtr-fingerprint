#include "IndigoCppFramework.h"
#include "molecule/molecule_fingerprint.h"
#include "indigo_internal.h"
#include "IndigoSubstructureMatcher.h"
#include "molecule/molecule_substructure_matcher.h"
#include "base_cpp/scanner.h"
#include "Profiling.h"

namespace {
    const std::shared_ptr<indigo_cpp::IndigoSession> globalIndigoSession = indigo_cpp::IndigoSession::create();

    std::unique_ptr<IndigoCppFramework::FingerprintT>
    tryBuildFingerprintFromMolecule(const IndigoCppFramework::MoleculeT &molecule) {
        INDIGO_BEGIN
                {
                    indigo::Molecule &mol = self.getObject(molecule.id()).getMolecule();
                    MoleculeFingerprintParameters fp_params(self.fp_params);
                    fp_params.ext = false;
                    fp_params.tau_qwords = 0;
//                    mol.aromatize(self.arom_options);
                    assert(mol.isAromatized());
                    indigo::MoleculeFingerprintBuilder fingerprintBuilder(mol, fp_params);
                    fingerprintBuilder.parseFingerprintType("sub", false);
                    fingerprintBuilder.process();
                    auto fpSize = fp_params.fingerprintSize();
                    assert(fpSize * CHAR_BIT == IndigoCppFramework::getFingerprintSize());

                    auto result = std::make_unique<indigo::Array<byte>>();
                    result->reserve(fpSize);
                    result->copy(fingerprintBuilder.get(), fpSize);
                    assert(result->size() * CHAR_BIT == IndigoCppFramework::getFingerprintSize());
                    return result;
                }
        INDIGO_END(nullptr);
    }

}

std::unique_ptr<IndigoCppFramework::MoleculeT> IndigoCppFramework::moleculeFromSmiles(const std::string &smiles) {

    auto mol = std::make_unique<MoleculeT>(globalIndigoSession->loadMolecule(smiles));
    mol->aromatize();
    return mol;
}

std::string IndigoCppFramework::moleculeToSmiles(const IndigoCppFramework::MoleculeT &molecule) {
    return molecule.smiles();
}

std::unique_ptr<IndigoCppFramework::QueryMoleculeT> IndigoCppFramework::queryMoleculeFromSmiles(const std::string &smiles) {
    auto mol = std::make_unique<QueryMoleculeT>(globalIndigoSession->loadQueryMolecule(smiles));
    mol->aromatize();
    return mol;
}

std::unique_ptr<IndigoCppFramework::FingerprintT>
IndigoCppFramework::fingerprintFromMolecule(const IndigoCppFramework::MoleculeT &molecule) {
    ProfileScope("IndigoCppFramework::fingerprintFromMolecule");
    auto res = tryBuildFingerprintFromMolecule(molecule);
    assert(res != nullptr);
    return res;
}

std::unique_ptr<IndigoCppFramework::StorageMoleculeT>
IndigoCppFramework::compressMolecule(const IndigoCppFramework::MoleculeT &molecule) {
    ProfileScope("IndigoCppFramework::compressMolecule");
    return std::make_unique<StorageMoleculeT>(molecule);
}

std::unique_ptr<IndigoCppFramework::MoleculeT>
IndigoCppFramework::decompressMolecule(const IndigoCppFramework::StorageMoleculeT &compressedMolecule) {
    ProfileScope("IndigoCppFramework::decompressMolecule");
    return std::make_unique<MoleculeT>(compressedMolecule);
}

bool IndigoCppFramework::isSubstructure(const IndigoCppFramework::QueryMoleculeT &queryMolecule,
                                        const IndigoCppFramework::MoleculeT &molecule) {
    ProfileScope("IndigoCppFramework::isSubstructure");
    Indigo& self = indigoGetInstance();
    indigo::Molecule &mol = self.getObject(molecule.id()).getMolecule();
    indigo::QueryMolecule &queryMol = self.getObject(queryMolecule.id()).getQueryMolecule();
    MoleculeSubstructureMatcher msm(mol);
    msm.setQuery(queryMol);
    bool res = msm.find();
    return res;
//    auto matcher = IndigoCppFramework::getGlobalIndigoSession()->substructureMatcher(molecule);
//    return matcher.match(queryMolecule);
}

size_t IndigoCppFramework::getFingerprintSize() {
    return 3072;
}

bool IndigoCppFramework::getFingerprintBit(const IndigoCppFramework::FingerprintT &fingerprint, size_t idx) {
    ProfileScope("IndigoCppFramework::getFingerprintBit");
    int byteIdx = int(idx / CHAR_BIT);
    int bitIndex = int(idx % CHAR_BIT);
    unsigned char byte = fingerprint.at(byteIdx);
    return bool((byte >> bitIndex) & 1u);
}

void IndigoCppFramework::setFingerprintBit(IndigoCppFramework::FingerprintT &fingerprint, size_t idx, bool val) {
    int byteIdx = int(idx / CHAR_BIT);
    int bitIndex = int(idx % CHAR_BIT);
    unsigned char &byte = fingerprint[byteIdx];

    if (val)
        byte |= (1 << bitIndex);
    else
        byte &= ~(1 << bitIndex);
}

bool IndigoCppFramework::isSubFingerprint(const IndigoCppFramework::FingerprintT &fingerprint1,
                                          const IndigoCppFramework::FingerprintT &fingerprint2) {
    assert(fingerprint1.size() == fingerprint2.size());
    ProfileScope("IndigoCppFramework::isSubFingerprint");
    for (int i = 0; i < fingerprint1.size(); i++) {
        unsigned char c1 = fingerprint1.at(i);
        unsigned char c2 = fingerprint2.at(i);
        if ((c1 & c2) != c1) {
            return false;
        }
    }
    return true;
}

std::shared_ptr<indigo_cpp::IndigoSession> IndigoCppFramework::getGlobalIndigoSession() {
    return globalIndigoSession;
}

IndigoCppFramework::FingerprintT IndigoCppFramework::getEmptyFingerprint() {
    FingerprintT fp;
    int byteSize = (int) getFingerprintSize() / CHAR_BIT;
    fp.reserve(byteSize);
    fp.resize(byteSize);
    fp.zerofill();
    return fp;
}


