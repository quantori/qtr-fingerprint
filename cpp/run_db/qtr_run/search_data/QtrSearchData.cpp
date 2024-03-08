#include "QtrSearchData.h"

#include "IndigoException.h"
#include "QtrRamSearchData.h"
#include "IndigoCFRamFilter.h"
#include "IndigoSmilesDriveFilter.h"
#include "PropertiesRamFilter.h"
#include "PropertiesDriveFilter.h"
#include "CompoundFilter.h"
#include "AlwaysTrueFilter.h"
#include "QtrDriveSearchData.h"


using namespace std;
using namespace indigo_cpp;

namespace qtr {
    namespace {
        unique_ptr <ByIdAnswerFilter>
        createRamFilter(const QtrRamSearchData &searchData, const std::string &querySmiles,
                        const PropertiesFilter::Bounds queryBounds) {
            auto indigoFilter = make_unique<IndigoCFRamFilter>(searchData.cfStorage, querySmiles);
            unique_ptr<ByIdAnswerFilter> filter = nullptr;
            if (searchData.propertiesTable != nullptr) {
                auto propertiesFilter = make_unique<PropertiesRamFilter>(searchData.propertiesTable,
                                                                         queryBounds);
                filter = make_unique<CompoundFilter<CIDType>>(std::move(propertiesFilter), std::move(indigoFilter));
            } else {
                filter = std::move(indigoFilter);
            }
            return std::move(filter);
        }

        unique_ptr <ByIdAnswerFilter>
        createDriveFilter(const QtrDriveSearchData &searchData, const std::string &querySmiles,
                          const PropertiesFilter::Bounds queryBounds) {
            auto querySmilesPtr = make_shared<string>(querySmiles);
            auto indigoFilter = make_unique<IndigoSmilesDriveFilter>(querySmilesPtr);
            auto propertiesFilter = make_unique<PropertiesDriveFilter>(queryBounds);
            auto filter = make_unique<CompoundFilter<CIDType>>(std::move(propertiesFilter), std::move(indigoFilter));
            return std::move(filter);
        }

        unique_ptr <ByIdAnswerFilter> createFilter(const QtrSearchData &searchData, const std::string &querySmiles,
                                                   const PropertiesFilter::Bounds queryBounds) {
            unique_ptr<ByIdAnswerFilter> filter;
            if (dynamic_cast<const QtrRamSearchData *>(&searchData) != nullptr) {
                const auto &ramSearchData = dynamic_cast<const QtrRamSearchData &>(searchData);
                filter = createRamFilter(ramSearchData, querySmiles, queryBounds);
            } else if (dynamic_cast<const QtrDriveSearchData *>(&searchData) != nullptr) {
                const auto &driveSearchData = dynamic_cast<const QtrDriveSearchData &>(searchData);
                filter = createDriveFilter(driveSearchData, querySmiles, queryBounds);
            } else {
                assert(false && "Undefined search data type");
            }
            return std::move(filter);
        }
    }

    QtrSearchData::QtrSearchData(shared_ptr<const BallTreeSearchEngine> ballTree,
                                 shared_ptr<const IdConverter> idConverter, size_t ansCount, size_t threadCount,
                                 double timeLimit, bool verificationStage) :
            SearchData(ansCount, threadCount, timeLimit, verificationStage),
            ballTree(std::move(ballTree)), idConverter(std::move(idConverter)) {}

    unique_ptr <QueryData<CIDType>>
    QtrSearchData::search(const SearchData::Query &query,
                          const PropertiesFilter::Bounds &queryBounds) {
        assert(query.smiles != nullptr);
        LOG(INFO) << "Start search: " << *query.smiles;
        Fingerprint fingerprint = getFingerprint(query);
        auto filter = getFilter(query, queryBounds);
        auto queryData = make_unique<BallTreeQueryData>(ansCount, timeLimit, fingerprint, std::move(filter));
        ballTree->search(*queryData, threadsCount);
        return std::move(queryData);
    }

    std::unique_ptr<ByIdAnswerFilter>
    QtrSearchData::getFilter(const SearchData::Query &query, const PropertiesFilter::Bounds &queryBounds) const {
        if (verificationStage) {
            return createFilter(*this, *query.smiles, queryBounds);
        } else
            return make_unique<AlwaysTrueFilter<CIDType>>();
    }

    Fingerprint QtrSearchData::getFingerprint(const SearchData::Query &query) {
        Fingerprint fingerprint;
        if (query.fingerprint != nullptr) {
            fingerprint = *query.fingerprint;
        } else {
            try {
                fingerprint = indigoFingerprintFromSmiles(*query.smiles);
            }
            catch (const std::exception &exception) {
                LOG(WARNING) << "Cannot build fingerprint: " << exception.what();
                return {};
            }
        }
        return fingerprint;
    }
} // qtr