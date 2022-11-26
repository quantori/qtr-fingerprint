#include "FingerprintTable.h"

#include "QtrIndigoFingerprint.h"

namespace qtr {

    template<>
    IndigoFingerprintTable::FingerprintTable(const std::string &sdfFile) {

        using namespace indigo_cpp;

        IndigoSessionPtr indigoSessionPtr = IndigoSession::create();
        IndigoSDFileIterator iterator = indigoSessionPtr->iterateSDFile(sdfFile);

        for (IndigoMoleculeSPtr &molecule: iterator) {
            molecule->aromatize();
            QtrIndigoFingerprint fingerprint(*molecule, "sub");
            this->emplace_back(fingerprint);
        }
    }

} // namespace qtr