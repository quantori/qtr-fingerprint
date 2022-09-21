#include "indigo.h"

#include "IndigoMolecule.h"
#include "IndigoSession.h"
#include "IndigoWriteBuffer.h"
#include "IndigoSDFileIterator.h"

#include "Fingerprint.h"

using namespace indigo_cpp;

namespace qtr {

    IndigoFingerprint cutFullIndigoFingerprint(const FullIndigoFingerprint &fullFingerprint) {
        IndigoFingerprint fingerprint;
        for (size_t i = 0; i < FullIndigoFingerprint::size(); i++) {
            if (i < FullIndigoFingerprintEmptyRange.first) {
                fingerprint[i] = fullFingerprint[i];
            } else if (i >= FullIndigoFingerprintEmptyRange.second) {
                fingerprint[i - FullIndigoFingerprintEmptyRange.second +
                            FullIndigoFingerprintEmptyRange.first] = fullFingerprint[i];
            }
        }
        return fingerprint;
    }

    IndigoFingerprint IndigoFingerprintFromSmiles(const std::string &smiles) {
        auto indigoSessionPtr = IndigoSession::create();
        auto mol = indigoSessionPtr->loadMolecule(smiles);
        mol.aromatize();
        int fingerprint = indigoFingerprint(mol.id(), "sub");
        FullIndigoFingerprint fullFingerprints(indigoToString(fingerprint));
        IndigoFingerprint cutFingerprint = cutFullIndigoFingerprint(fullFingerprints);
        return cutFingerprint;
    }
}