#include "Fingerprint.h"

namespace qtr {

    Fingerprint::Fingerprint(size_t size) : Bitset<>(size) {}

    Fingerprint::Fingerprint(const std::string &s) : Bitset<>(s.size() * BIT_IN_BYTE / 2) {
        assert(s.size() % 2 == 0);
        assert(s.size() == sizeInBytes() * 2);
        for (size_t i = 0; i < s.size(); ++i) {
            addSymbol(i, 4, (size_t) chexToInt(s[i]));
        }
    }

    Fingerprint::Fingerprint(const Array<byte> &arr) : Bitset<>(arr.size() * BIT_IN_BYTE) {
        assert(arr.size() == sizeInBytes());
        for (size_t i = 0; i < arr.size(); i++) {
            addSymbol(i, 8, arr[i]);
        }
    }

    void Fingerprint::setBytes(const std::vector<std::byte> &bytes) {
        assert(bytes.size() == sizeInBytes());
        for (size_t i = 0; i < bytes.size(); i++)
            for (size_t j = 0; j < CHAR_BIT; j++)
                this->operator[](i * CHAR_BIT + j) = bool((bytes[i] >> j) & std::byte(1));
    }

    std::vector<std::byte> Fingerprint::getBytes() const {
        std::vector<std::byte> result(sizeInBytes(), std::byte(0));
        for (size_t i = 0; i < result.size(); i++)
            for (size_t j = 0; j < CHAR_BIT; j++)
                if (this->operator[](i * CHAR_BIT + j))
                    result[i] |= (std::byte(1) << j);

        return result;
    }

    void Fingerprint::addSymbol(size_t blockIndex, size_t blockSize, size_t symbol) {
        for (size_t i = 0; i < blockSize; i++) {
            (*this)[blockIndex * blockSize + blockSize - i - 1] = symbol & (1u << i);
        }
    }

    Fingerprint::Fingerprint(const QtrIndigoFingerprint &f) : Bitset<>(IndigoFingerprintSize) {
        setBytes(f.data());
    }

//    Fingerprint::Fingerprint(std::unique_ptr<ExplicitBitVect> f) {
    // TODO

//    }

    Fingerprint cutFullIndigoFingerprint(const Fingerprint &fullFingerprint) {
        assert(fullFingerprint.size() == FullIndigoFingerprintSize);
        Fingerprint fingerprint(IndigoFingerprintSize);
        for (size_t i = 0; i < FullIndigoFingerprintSize; i++) {
            if (i < FullIndigoFingerprintEmptyRange.first) {
                fingerprint[i] = fullFingerprint[i];
            } else if (i >= FullIndigoFingerprintEmptyRange.second) {
                fingerprint[i - FullIndigoFingerprintEmptyRange.second +
                            FullIndigoFingerprintEmptyRange.first] = fullFingerprint[i];
            } else {
                assert(!fullFingerprint[i]);
            }
        }
        return fingerprint;
    }

    const static Indigo indigoInstance;

    Fingerprint indigoFingerprintFromSmiles(const std::string &smiles) {
        indigo::BufferScanner scanner(smiles.c_str(), smiles.size(), false);
        indigo::SmilesLoader loader(scanner);
        indigo::Molecule molecule;
        loader.loadMolecule(molecule);
        molecule.aromatize(indigo::AromaticityOptions());
        bingo::IndexMolecule indexMolecule(molecule, indigo::AromaticityOptions());
        indigo::Array<byte> subFingerprint;
        indigo::MoleculeFingerprintBuilder fingerprintBuilder(molecule, indigoInstance.fp_params);
        fingerprintBuilder.parseFingerprintType("sub", false);
        fingerprintBuilder.process();
        subFingerprint.copy(fingerprintBuilder.get(), indigoInstance.fp_params.fingerprintSize());
//        indexMolecule.buildFingerprint(indigoInstance.fp_params, &subFingerprint, nullptr);
        Fingerprint fullIndigoFingerprint(subFingerprint);
        Fingerprint cutFingerprint = cutFullIndigoFingerprint(fullIndigoFingerprint);
        return cutFingerprint;
    }

    Fingerprint rdkitFingerprintFromSmiles(const std::string &smiles) {
        // TODO
        //        std::shared_ptr<RDKit::ROMol> mol3( RDKit::SmilesToMol( "Cc1cccc" ) );
        //        std::unique_ptr<ExplicitBitVect> mfp(PatternFingerprintMol(*mol));
        //        Cast ExplicitBitVect to Fingerprint
        return Fingerprint();
    }
};