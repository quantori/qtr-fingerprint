#include "BingoNoSQLSearchData.h"

#include "QueryData.h"

#include <utility>

#include "IndigoException.h"
#include "BingoNoSQL.h"


using namespace std;
using namespace indigo_cpp;

namespace qtr {
    namespace {
        void searchInRDKitDB(QueryData<CIDType> &queryData, uint64_t ansCount, BingoMolecule &db,
                             const IndigoQueryMolecule &query) {
            vector<CIDType> result;
            auto subMatcher = db.searchSub(query, "");
            for (auto &mol: subMatcher) {
                if (result.size() >= ansCount || queryData.checkShouldStop() || queryData.checkTimeOut())
                    break;
                result.push_back(mol.getId());
            }
            if (queryData.checkTimeOut()) {
                LOG(INFO) << "Search stopped due to timeout";
            }
            queryData.addAnswers(result);
            queryData.tagFinishTask();
        }
    }

    BingoNoSQLSearchData::BingoNoSQLSearchData(const std::filesystem::path &dbDataDir, size_t ansCount,
                                               size_t threadsCount, double timeLimit,
                                               bool verificationStage) :
            SearchData(ansCount, threadsCount, timeLimit, verificationStage),
            db(BingoMolecule::loadDatabaseFile(IndigoSession::create(), dbDataDir, "")) {
        if (!verificationStage) {
            throw std::runtime_error("VerificationStage=false is not supported for BingoNoSQL");
        }
    }

    unique_ptr<QueryData<CIDType>>
    BingoNoSQLSearchData::search(const SearchData::Query &query, const PropertiesFilter::Bounds &) {
        assert(query.smiles != nullptr);
        assert(query.fingerprint == nullptr);
        LOG(INFO) << "Start search: " << *query.smiles;
        unique_ptr<IndigoQueryMolecule> molecule = nullptr;
        try {
            molecule = make_unique<IndigoQueryMolecule>(db.session->loadQueryMolecule(*query.smiles));
            molecule->aromatize();
        }
        catch (const IndigoException &e) {
            LOG(ERROR) << "Cannot parse smiles: " << *query.smiles << " (" << e.what() << ")";
            return nullptr;
        }

        auto queryData = make_unique<QueryData<CIDType>>(ansCount, timeLimit, make_unique<AlwaysTrueFilter<CIDType>>());
        queryData->addTask(std::async(launch::async, searchInRDKitDB, ref(*queryData), ansCount, std::ref(db),
                                      std::move(*molecule)));
        return std::move(queryData);
    }

    BingoNoSQLSearchData::~BingoNoSQLSearchData() {
        db.close();
    }
} // qtr