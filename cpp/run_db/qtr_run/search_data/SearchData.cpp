#include "SearchData.h"

namespace qtr {
    SearchData::SearchData(size_t ansCount, size_t threadCount, double timeLimit, bool verificationStage) :
            ansCount(ansCount), threadsCount(threadCount), timeLimit(timeLimit), verificationStage(verificationStage) {}

    BaseLibrary SearchData::getBaseLibrary() const {
        return BaseLibrary::BadOption;
    }

    Fingerprint SearchData::Query::getFingerprint() const {
        if (fingerprint != nullptr) {
            return *fingerprint;
        } else {
            try {
                if (baseLibrary == BaseLibrary::RDKit) {
                    return rdkitFingerprintFromSmiles(*smiles);
                } else if (baseLibrary == BaseLibrary::Indigo) {
                    return indigoFingerprintFromSmiles(*smiles);
                } else {
                    throw std::runtime_error("Invalid baseLibrary value");
                }
            }
            catch (const std::exception &exception) {
                LOG_ERROR_AND_EXIT(std::string("Cannot build fingerprint for ") + *smiles + " : " + exception.what());
            }
        }
    }
} // qtr
