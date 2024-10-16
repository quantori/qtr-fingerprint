#include "BallTreeSearchData.h"

#include "IndigoException.h"
#include "QtrRamSearchData.h"
#include "IndigoCFRamFilter.h"
#include "IndigoSmilesDriveFilter.h"
#include "PropertiesRamFilter.h"
#include "PropertiesDriveFilter.h"
#include "CompoundFilter.h"
#include "AlwaysTrueFilter.h"
#include "QtrDriveSearchData.h"
#include "RDKitSubstructFilter.h"


using namespace std;
using namespace indigo_cpp;

namespace qtr {
    namespace {
        unique_ptr <ByIdAnswerFilter>
        getSubstructFilter(const QtrRamSearchData &searchData, const string &querySmiles) {
            unique_ptr<ByIdAnswerFilter> substructFilter;
            assert(searchData.cfStorage == nullptr ^ searchData.molHolder == nullptr);
            if (searchData.molHolder != nullptr) {
                assert(searchData.cfStorage == nullptr);
                substructFilter = make_unique<RDKitSubstructFilter>(searchData.molHolder, querySmiles);
            } else if (searchData.cfStorage != nullptr) {
                assert(searchData.molHolder == nullptr);
                substructFilter = make_unique<IndigoCFRamFilter>(searchData.cfStorage, querySmiles);
            } else {
                LOG_ERROR_AND_EXIT("Cannot create Ram Filter because no molecules storage specified");
            }
            return substructFilter;
        }

        unique_ptr <ByIdAnswerFilter>
        addPropertiesToFilter(unique_ptr <ByIdAnswerFilter> filter, const QtrRamSearchData &searchData,
                              const PropertiesFilter::Bounds queryBounds) {
            if (searchData.propertiesTable != nullptr) {
                auto propertiesFilter = make_unique<PropertiesRamFilter>(searchData.propertiesTable,
                                                                         queryBounds);
                return make_unique<CompoundFilter<CIDType>>(std::move(propertiesFilter), std::move(filter));
            } else {
                return filter;
            }
        }

        unique_ptr <ByIdAnswerFilter>
        createRamFilter(const QtrRamSearchData &searchData, const std::string &querySmiles,
                        const PropertiesFilter::Bounds queryBounds) {
            auto substructFilter = getSubstructFilter(searchData, querySmiles);
            return addPropertiesToFilter(std::move(substructFilter), searchData, queryBounds);
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

        unique_ptr <ByIdAnswerFilter> createFilter(const BallTreeSearchData &searchData, const std::string &querySmiles,
                                                   const PropertiesFilter::Bounds queryBounds) {
            unique_ptr<ByIdAnswerFilter> filter;
            if (dynamic_cast<const QtrRamSearchData *>(&searchData) != nullptr) {
                const auto &ramSearchData = dynamic_cast<const QtrRamSearchData &>(searchData);
                filter = createRamFilter(ramSearchData, querySmiles, queryBounds);
            } else if (dynamic_cast<const QtrDriveSearchData *>(&searchData) != nullptr) {
                const auto &driveSearchData = dynamic_cast<const QtrDriveSearchData &>(searchData);
                filter = createDriveFilter(driveSearchData, querySmiles, queryBounds);
            } else {
                LOG_ERROR_AND_EXIT("Undefined search data type");
            }
            return std::move(filter);
        }
    }

    BallTreeSearchData::BallTreeSearchData(shared_ptr<const BallTreeSearchEngine> ballTree,
                                           shared_ptr<const IdConverter> idConverter, size_t ansCount,
                                           size_t threadCount,
                                           double timeLimit, bool verificationStage) :
            SearchData(ansCount, threadCount, timeLimit, verificationStage),
            ballTree(std::move(ballTree)), idConverter(std::move(idConverter)) {}

    unique_ptr <QueryData<CIDType>>
    BallTreeSearchData::search(const SearchData::Query &query,
                               const PropertiesFilter::Bounds &queryBounds) {
        assert(query.smiles != nullptr);
        LOG(INFO) << "Start search: " << *query.smiles;
        Fingerprint fingerprint = query.getFingerprint();
        if (fingerprint.size() == 0) {
            return nullptr;
        }
        auto filter = getFilter(query, queryBounds);
        auto queryData = make_unique<QueryDataWithFingerprint>(ansCount, timeLimit, fingerprint, std::move(filter));
        ballTree->search(*queryData, threadsCount);
        return std::move(queryData);
    }

    std::unique_ptr<ByIdAnswerFilter>
    BallTreeSearchData::getFilter(const SearchData::Query &query, const PropertiesFilter::Bounds &queryBounds) const {
        if (verificationStage) {
            return createFilter(*this, *query.smiles, queryBounds);
        } else
            return make_unique<AlwaysTrueFilter<CIDType>>();
    }
} // qtr