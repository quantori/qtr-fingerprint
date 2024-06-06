#include "QtrRamSearchData.h"

#include <utility>

namespace qtr {
    QtrRamSearchData::QtrRamSearchData(std::shared_ptr<const BallTreeSearchEngine> ballTree,
                                       std::shared_ptr<const IdConverter> idConverter, size_t ansCount,
                                       size_t threadCount, double timeLimit, std::shared_ptr<CFStorage> cfStorage,
                                       std::shared_ptr<RDKit::MolHolderBase> molHolder,
                                       std::shared_ptr<const std::vector<PropertiesFilter::Properties>> propertiesTable,
                                       bool verificationStage)
            : BallTreeSearchData(std::move(ballTree), std::move(idConverter), ansCount, threadCount, timeLimit,
                                 verificationStage), cfStorage(std::move(cfStorage)), molHolder(std::move(molHolder)),
              propertiesTable(std::move(propertiesTable)) {
    }

    BaseLibrary QtrRamSearchData::getBaseLibrary() const {
        return cfStorage != nullptr ? BaseLibrary::Indigo :
               molHolder != nullptr ? BaseLibrary::RDKit : BaseLibrary::BadOption;
    }
} // qtr