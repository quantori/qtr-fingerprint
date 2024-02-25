#include "FingerprintTable.h"

#include "QtrIndigoFingerprint.h"
#include "IndigoIterator.h"

namespace qtr {

    FingerprintTable::FingerprintTable(const std::string &sdfFile) {

        using namespace indigo_cpp;

        IndigoSessionPtr indigoSessionPtr = IndigoSession::create();
        auto iterator = indigoSessionPtr->iterateSDFile(sdfFile);

        for (IndigoMoleculeSPtr &molecule: iterator) {
            molecule->aromatize();
            QtrIndigoFingerprint fingerprint(*molecule, "sub");
            this->emplace_back(fingerprint);
        }
    }

} // namespace qtr