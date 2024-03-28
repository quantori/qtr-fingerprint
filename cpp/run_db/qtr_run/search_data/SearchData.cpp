#include "SearchData.h"

namespace qtr {
    SearchData::SearchData(size_t ansCount, size_t threadCount, double timeLimit, bool verificationStage) :
            ansCount(ansCount), threadsCount(threadCount), timeLimit(timeLimit), verificationStage(verificationStage) {}

    Fingerprint SearchData::Query::getFingerprint() const {
        Fingerprint res;
        if (fingerprint != nullptr) {
            res = *fingerprint;
        } else {
            try {
                res = indigoFingerprintFromSmiles(*smiles);
            }
            catch (const std::exception &exception) {
                LOG(WARNING) << "Cannot build fingerprint: " << exception.what();
                return {};
            }
        }
        return res;
    }
} // qtr
