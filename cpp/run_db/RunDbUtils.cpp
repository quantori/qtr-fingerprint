#include "RunDbUtils.h"

#include "string_table_io/StringTableReader.h"
#include "BallTreeQueryData.h"
#include "IndigoQueryMolecule.h"
#include "IndigoSubstructureMatcher.h"
#include "IndigoRamSmilesFilter.h"
#include "IndigoRamMoleculeFilter.h"
#include "IndigoDriveFilter.h"
#include "PropertiesRamFilter.h"
#include "PropertiesDriveFilter.h"
#include "CompoundFilter.h"
#include "search_data/RamSmilesSearchData.h"
#include "search_data/DriveSearchData.h"
#include "search_data/RamMoleculesSearchData.h"


using namespace std;

namespace qtr {

    namespace {
        unique_ptr <AnswerFilter>
        createRamSmilesFilter(const RamSmilesSearchData &searchData, const string &querySmiles,
                              const PropertiesFilter::Bounds queryBounds) {
            auto querySmilesPtr = make_shared<string>(querySmiles);
            auto indigoFilter = make_unique<IndigoRamSmilesFilter>(searchData.smilesTable, querySmilesPtr);
            auto propertiesFilter = make_unique<PropertiesRamFilter>(searchData.propertiesTable,
                                                                     queryBounds);
            auto filter = make_unique<CompoundFilter>(std::move(propertiesFilter), std::move(indigoFilter));
            return std::move(filter);
        }

        unique_ptr <AnswerFilter>
        createRamMoleculesFilter(const RamMoleculesSearchData &searchData, const std::string &querySmiles,
                                 const PropertiesFilter::Bounds queryBounds) {
            auto querySmilesPtr = make_shared<string>(querySmiles);
            auto indigoFilter = make_unique<IndigoRamMoleculeFilter>(searchData.moleculesTable, querySmilesPtr);
            auto propertiesFilter = make_unique<PropertiesRamFilter>(searchData.propertiesTable,
                                                                     queryBounds);
            auto filter = make_unique<CompoundFilter>(std::move(propertiesFilter), std::move(indigoFilter));
            return std::move(filter);
        }

        unique_ptr <AnswerFilter> createDriveFilter(const DriveSearchData &searchData, const string &querySmiles,
                                                    const PropertiesFilter::Bounds queryBounds) {
            auto querySmilesPtr = make_shared<string>(querySmiles);
            auto indigoFilter = make_unique<IndigoDriveFilter>(querySmilesPtr);
            auto propertiesFilter = make_unique<PropertiesDriveFilter>(queryBounds);
            auto filter = make_unique<CompoundFilter>(std::move(propertiesFilter), std::move(indigoFilter));
            return std::move(filter);
        }

        unique_ptr <AnswerFilter> createFilter(const SearchData &searchData, const std::string &querySmiles,
                                               const PropertiesFilter::Bounds queryBounds) {
            unique_ptr<AnswerFilter> filter;
            if (searchData.getDbType() == DbType::InRamSmiles) {
                const auto &ramSearchData = dynamic_cast<const RamSmilesSearchData &>(searchData);
                filter = createRamSmilesFilter(ramSearchData, querySmiles, queryBounds);
            } else if (searchData.getDbType() == DbType::OnDrive) {
                const auto &driveSearchData = dynamic_cast<const DriveSearchData &>(searchData);
                filter = createDriveFilter(driveSearchData, querySmiles, queryBounds);
            } else if (searchData.getDbType() == DbType::InRamMolecules) {
                const auto &ramMoleculesSearchData = dynamic_cast<const RamMoleculesSearchData &>(searchData);
                filter = createRamMoleculesFilter(ramMoleculesSearchData, querySmiles, queryBounds);
            } else {
                assert(false && "Undefined search data type");
            }
            return std::move(filter);
        }
    }

    pair<bool, unique_ptr<BallTreeQueryData>>
    runSearch(const SearchData &searchData, const string &querySmiles, const PropertiesFilter::Bounds &queryBounds) {
        LOG(INFO) << "Start search: " << querySmiles;
        IndigoFingerprint fingerprint;
        try {
            fingerprint = qtr::indigoFingerprintFromSmiles(querySmiles);
        }
        catch (exception &exception) {
            LOG(WARNING) << "Skip query:" << exception.what();
            return {true, unique_ptr<BallTreeQueryData>(nullptr)};
        }
        auto filter = createFilter(searchData, querySmiles, queryBounds);
        auto queryData = make_unique<BallTreeQueryData>(searchData.ansCount, fingerprint, std::move(filter));
        searchData.ballTree->search(*queryData, searchData.threadsCount);
        return {false, std::move(queryData)};
    }

} // namespace qtr

