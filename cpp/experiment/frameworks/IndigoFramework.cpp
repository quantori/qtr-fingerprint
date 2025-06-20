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

MoleculeFingerprintParameters IndigoFramework::getFingerprintParams() const {
    auto &self = indigoGetInstance();
    MoleculeFingerprintParameters params(self.fp_params);
    params.ext = true;
    params.tau_qwords *= 0;
    params.sim_qwords *= 0;
    params.ord_qwords *= _fingerprintRatio;
    params.any_qwords *= 0;
    return params;
}

std::unique_ptr<IndigoFramework::FingerprintT>
IndigoFramework::tryBuildFingerprintFromMolecule(const IndigoFramework::MoleculeT &molecule) const {
    INDIGO_BEGIN
            {
                auto &mol = const_cast<Molecule &>(molecule);
                auto fpParams = getFingerprintParams();
//                    mol.aromatize(self.arom_options);
                assert(mol.isAromatized());
                indigo::MoleculeFingerprintBuilder fingerprintBuilder(mol, fpParams);
                fingerprintBuilder.parseFingerprintType("sub", false);
                fingerprintBuilder.process();
                auto fpSizeBytes = fpParams.fingerprintSize();

            int sizeRate = sizeof(IndigoFramework::FingerprintInnerT) / sizeof(byte);
            auto fpSize = (fpSizeBytes + sizeRate - 1) / sizeRate;
            assert(fpSize % sizeRate == 0);

                auto result = std::make_unique<indigo::Array<IndigoFramework::FingerprintInnerT>>();
                result->reserve(fpSize);
                result->copy(reinterpret_cast<const IndigoFramework::FingerprintInnerT *>(fingerprintBuilder.get()),
                                fpSize);
                return result;
            }
    INDIGO_END(nullptr);
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
}

std::unique_ptr<IndigoFramework::FingerprintT>
IndigoFramework::fingerprintFromMolecule(const IndigoFramework::MoleculeT &molecule) const {
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
//    Indigo &self = indigoGetInstance();
    auto &mutableMol = const_cast<IndigoFramework::MoleculeT &>(molecule);
    auto &mutableQueryMol = const_cast<IndigoFramework::QueryMoleculeT &>(queryMolecule);
    MoleculeSubstructureMatcher msm(mutableMol);
    msm.setQuery(mutableQueryMol);
    bool res = msm.find();
    return res;
}

size_t IndigoFramework::getFingerprintSize() const {
    return getFingerprintParams().fingerprintSize() * 8;
//                     ext  tau  sim  ord   any
//    return 3072;  //  -    -    +    +     +   (BingoNOSQL)
//    return 3096;  //  +    -    +    +     +
//    return 2584;  //  +    -    -    +     +
//    return 2560;  //  -    -    -    +     +
//    return 3224;  //  +    +    -    +     +
//    return 4184;  //  +    -    -    2x    +
//    return 12824; //  +    -    -    5x    5x
//    return 8024;  //  +    -    -    5x    -
//    return 12824; //  +    -    -    8x    -
}

bool IndigoFramework::getFingerprintBit(const IndigoFramework::FingerprintT &fingerprint, size_t idx) {
//    ProfileScope("IndigoFramework::getFingerprintBit");
    constexpr int itemBitSize = sizeof(IndigoFramework::FingerprintInnerT) * 8;
    int byteIdx = int(idx / itemBitSize);
    int bitIndex = int(idx % itemBitSize);
    unsigned char byte = fingerprint.at(byteIdx);
    return bool((byte >> bitIndex) & 1u);
}

void IndigoFramework::setFingerprintBit(IndigoFramework::FingerprintT &fingerprint, size_t idx, bool val) {
    constexpr int itemBitSize = sizeof(IndigoFramework::FingerprintInnerT) * 8;
    int byteIdx = int(idx / itemBitSize);
    int bitIndex = int(idx % itemBitSize);
    auto &item = fingerprint[byteIdx];

    if (val)
        item |= (1 << bitIndex);
    else
        item &= ~(1 << bitIndex);
}

bool IndigoFramework::isSubFingerprint(const IndigoFramework::QueryFingerprintT &fingerprint1,
                                       const IndigoFramework::FingerprintT &fingerprint2) {
    ProfileScope("IndigoFramework::isSubFingerprint");
    return fingerprint1.isSubFingerprint(fingerprint2);
}

IndigoFramework::FingerprintT IndigoFramework::getEmptyFingerprint() const {
    FingerprintT fp;
    constexpr int itemBitSize = sizeof(IndigoFramework::FingerprintInnerT) * 8;
    int fpSize = (int) (getFingerprintSize() + itemBitSize - 1) / itemBitSize;
    fp.reserve(fpSize);
    fp.resize(fpSize);
    fp.zerofill();
    return fp;
}

std::unique_ptr<IndigoFramework::QueryFingerprintT>
IndigoFramework::queryFingerprintFromFingerprint(const IndigoFramework::FingerprintT &fingerprint) {
    return std::make_unique<IndigoFramework::QueryFingerprintT>(fingerprint);
}

void IndigoFramework::init(const Config &config) {
    _fingerprintRatio = config.getDouble("fpRatio", 1.0);
}

IndigoFramework &IndigoFramework::getInstance() {
    static IndigoFramework instance;
    return instance;
}
