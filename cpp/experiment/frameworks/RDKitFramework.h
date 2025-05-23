#pragma once

#include <string>

#include "GraphMol/ROMol.h"
#include "GraphMol/Fingerprints/Fingerprints.h"
#include "GraphMol/Substruct/SubstructMatch.h"
#include "DataStructs/ExplicitBitVect.h"

#include "frameworks/FrameworkInterface.h"
#include "Bitset.h"

template<typename T>
class RDKitQueryFingerprint {
public:
    explicit RDKitQueryFingerprint(const Bitset<T> &fingerprint) {
        const auto *data = fingerprint.data();
        for (size_t i = 0; i < fingerprint.dataVectorSize(); i++) {
            if (data[i] == 0) {
                continue;
            }
            _items.emplace_back(i, data[i]);
        }
    }

    [[nodiscard]] bool isSubFingerprint(const Bitset<T> &fingerprint) const {
        const auto *otherData = fingerprint.data();
        for (auto& [i, data]: _items) {
            if ((data & otherData[i]) != data) {
                return false;
            }
        }
        return true;
    }

private:
    std::vector<std::pair<size_t, T>> _items;
};

class RDKitFramework {
public:
    using MoleculeT = RDKit::ROMol;
    using QueryMoleculeT = RDKit::ROMol;
    using StorageMoleculeT = std::string;

    using FingerprintT = Bitset<unsigned long long>;
    using QueryFingerprintT = RDKitQueryFingerprint<unsigned long long>;

    static std::unique_ptr<MoleculeT> moleculeFromSmiles(const std::string &smiles);

    static std::string moleculeToSmiles(const MoleculeT &molecule);

    static std::unique_ptr<QueryMoleculeT> queryMoleculeFromSmiles(const std::string &smiles);

    static std::unique_ptr<FingerprintT> fingerprintFromMolecule(const MoleculeT &molecule);

    static std::unique_ptr<QueryFingerprintT> queryFingerprintFromFingerprint(const FingerprintT &fingerprint);

    static std::unique_ptr<StorageMoleculeT> compressMolecule(const MoleculeT &molecule);

    static std::unique_ptr<MoleculeT> decompressMolecule(const StorageMoleculeT &compressedMolecule);

    static bool isSubstructure(const QueryMoleculeT &queryMolecule, const MoleculeT &molecule);

    static RDKit::SubstructMatchParameters getSubstructMatchParameters();

    static bool getFingerprintBit(const FingerprintT &fingerprint, size_t idx);

    static size_t getFingerprintSize();

    static void setFingerprintBit(FingerprintT &fingerprint, size_t idx, bool val);

    static bool isSubFingerprint(const QueryFingerprintT &fingerprint1, const FingerprintT &fingerprint2);

    static FingerprintT getEmptyFingerprint();
};

static_assert(FrameworkInterface<RDKitFramework>, "RDKitFramework must satisfy FrameworkInterface");
