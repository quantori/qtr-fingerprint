#include "BingoNoSQLSearchData.h"

#include "QueryData.h"

#include <utility>

#include "IndigoException.h"
#include "BingoNoSQL.h"


using namespace std;
using namespace indigo_cpp;

namespace qtr {
    namespace {
        void searchInBingoDB(QueryData<CIDType> &queryData, uint64_t ansCount, BingoMolecule &db,
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
                                               size_t threadsCount, double timeLimit) :
            SearchData(ansCount, threadsCount, timeLimit),
            db(BingoMolecule::loadDatabaseFile(IndigoSession::create(), dbDataDir, "")) {}

    unique_ptr <QueryData<CIDType>>
    BingoNoSQLSearchData::search(const string &querySmiles, const PropertiesFilter::Bounds &) {
        LOG(INFO) << "Start search: " << querySmiles;
        unique_ptr<IndigoQueryMolecule> molecule = nullptr;
        try {
            molecule = make_unique<IndigoQueryMolecule>(db.session->loadQueryMolecule(querySmiles));
            molecule->aromatize();
        }
        catch (const IndigoException &e) {
            LOG(ERROR) << "Cannot parse smiles: " << querySmiles << " (" << e.what() << ")";
            return nullptr;
        }

        auto queryData = make_unique<QueryData<CIDType>>(ansCount, timeLimit, make_unique<AlwaysTrueFilter<CIDType>>());
        queryData->addTask(std::async(launch::async, searchInBingoDB, ref(*queryData), ansCount, std::ref(db),
                                      std::move(*molecule)));
        return std::move(queryData);
    }

    BingoNoSQLSearchData::~BingoNoSQLSearchData() {
        db.close();
    }
} // qtr