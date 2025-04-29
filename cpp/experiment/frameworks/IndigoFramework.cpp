#include "IndigoFramework.h"
#include "molecule/molecule_fingerprint.h"
#include "indigo_internal.h"
#include "IndigoSubstructureMatcher.h"
#include "molecule/molecule_substructure_matcher.h"
#include "base_cpp/scanner.h"
#include "Profiling.h"
#include "molecule/smiles_loader.h"
#include "molecule/cmf_loader.h"
#include "molecule/query_molecule.h"
#include "molecule/molecule.h"

namespace {
    std::unique_ptr<IndigoFramework::FingerprintT>
    tryBuildFingerprintFromMolecule(const IndigoFramework::MoleculeT &molecule) {
        INDIGO_BEGIN
                {
//                    indigo::Molecule &mol = self.getObject(molecule).getMolecule();
                    Molecule &mol = const_cast<Molecule &>(molecule);
                    MoleculeFingerprintParameters fp_params(self.fp_params);
                    fp_params.ext = false;
                    fp_params.tau_qwords = 0;
//                    mol.aromatize(self.arom_options);
                    assert(mol.isAromatized());
                    indigo::MoleculeFingerprintBuilder fingerprintBuilder(mol, fp_params);
                    fingerprintBuilder.parseFingerprintType("sub", false);
                    fingerprintBuilder.process();
                    auto fpSize = fp_params.fingerprintSize();
                    assert(fpSize * CHAR_BIT == IndigoFramework::getFingerprintSize());

                    auto result = std::make_unique<indigo::Array<byte>>();
                    result->reserve(fpSize);
                    result->copy(fingerprintBuilder.get(), fpSize);
                    assert(result->size() * CHAR_BIT == IndigoFramework::getFingerprintSize());
                    return result;
                }
        INDIGO_END(nullptr);
    }
}

std::unique_ptr<IndigoFramework::MoleculeT> IndigoFramework::moleculeFromSmiles(const std::string &smiles) {
    Indigo &self = indigoGetInstance();
    BufferScanner scanner(smiles.c_str());
    SmilesLoader loader(scanner);
    loader.stereochemistry_options = self.stereochemistry_options;
    loader.ignore_bad_valence = self.ignore_bad_valence;
    loader.ignore_no_chiral_flag = self.ignore_no_chiral_flag;

    auto molecule = std::make_unique<Molecule>();
    loader.loadMolecule(*molecule);
    molecule->aromatize(self.arom_options);
    return molecule;
}

std::string IndigoFramework::moleculeToSmiles(const IndigoFramework::MoleculeT &molecule) {
    std::string smiles;
    StringOutput output(smiles);
    SmilesSaver saver(output);
    Indigo &self = indigoGetInstance();
//    Molecule &mol = self.getObject(molecule).getMolecule();
    Molecule &mol = const_cast<Molecule &>(molecule);
    saver.saveMolecule(mol);
    return smiles;
}

std::unique_ptr<IndigoFramework::QueryMoleculeT> IndigoFramework::queryMoleculeFromSmiles(const std::string &smiles) {
    Indigo &self = indigoGetInstance();
    BufferScanner scanner(smiles.c_str());
    SmilesLoader loader(scanner);
    loader.stereochemistry_options = self.stereochemistry_options;
    loader.ignore_bad_valence = self.ignore_bad_valence;
    loader.ignore_no_chiral_flag = self.ignore_no_chiral_flag;

    auto molecule = std::make_unique<QueryMolecule>();
    loader.loadQueryMolecule(*molecule);
    molecule->aromatize(self.arom_options);
    return molecule;
//    int id = self.addObject(std::move(queryMolPtr));
//    return std::make_unique<IndigoFramework::QueryMoleculeT>(id);
}

std::unique_ptr<IndigoFramework::FingerprintT>
IndigoFramework::fingerprintFromMolecule(const IndigoFramework::MoleculeT &molecule) {
    ProfileScope("IndigoFramework::fingerprintFromMolecule");
    auto res = tryBuildFingerprintFromMolecule(molecule);
    assert(res != nullptr);
    return res;
}

std::unique_ptr<IndigoFramework::StorageMoleculeT>
IndigoFramework::compressMolecule(const IndigoFramework::MoleculeT &molecule) {
    ProfileScope("IndigoFramework::compressMolecule");
    Indigo &self = indigoGetInstance();
    auto result = std::make_unique<IndigoFramework::StorageMoleculeT>();
    StringOutput output(*result);
    CmfSaver saver(output);
//    Molecule &mol = self.getObject(molecule);
//    Molecule &mol = const_cast<Molecule &>(molecule);
    Molecule &mutableMol = const_cast<Molecule &>(molecule);
    saver.saveMolecule(mutableMol);
    return result;
}

std::unique_ptr<IndigoFramework::MoleculeT>
IndigoFramework::decompressMolecule(const IndigoFramework::StorageMoleculeT &compressedMolecule) {
    ProfileScope("IndigoFramework::decompressMolecule");
    BufferScanner scanner(compressedMolecule.c_str(), compressedMolecule.size(), false);
    CmfLoader cmfLoader(scanner);

    auto molecule = std::make_unique<IndigoFramework::MoleculeT>();
    cmfLoader.loadMolecule(*molecule);
    return molecule;
}

bool IndigoFramework::isSubstructure(const IndigoFramework::QueryMoleculeT &queryMolecule,
                                     const IndigoFramework::MoleculeT &molecule) {
    ProfileScope("IndigoFramework::isSubstructure");
    Indigo &self = indigoGetInstance();
    auto &mutableMol = const_cast<IndigoFramework::MoleculeT &>(molecule);
    auto &mutableQueryMol = const_cast<IndigoFramework::QueryMoleculeT &>(queryMolecule);
    MoleculeSubstructureMatcher msm(mutableMol);
    msm.setQuery(mutableQueryMol);
    bool res = msm.find();
    return res;
}

size_t IndigoFramework::getFingerprintSize() {
    return 3072;
}

bool IndigoFramework::getFingerprintBit(const IndigoFramework::FingerprintT &fingerprint, size_t idx) {
    ProfileScope("IndigoFramework::getFingerprintBit");
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

IndigoFramework::FingerprintT IndigoFramework::getEmptyFingerprint() {
    FingerprintT fp;
    int byteSize = (int) getFingerprintSize() / CHAR_BIT;
    fp.reserve(byteSize);
    fp.resize(byteSize);
    fp.zerofill();
    return fp;
}
