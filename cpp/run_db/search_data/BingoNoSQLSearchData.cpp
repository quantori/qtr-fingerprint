#include "BingoNoSQLSearchData.h"

#include "QueryData.h"

#include <utility>

#include "bingo-nosql.h"
#include "IndigoQueryMolecule.h"
#include "IndigoException.h"

using namespace std;
using namespace indigo_cpp;

namespace qtr {
    namespace {
        void searchInBingoDB(QueryData<CIDType> &queryData, uint64_t ansCount, const IndigoSessionPtr &indigoSessionPtr,
                             int db, const IndigoQueryMolecule &query) {
            queryData.tagStartTask();
            vector<CIDType> result;
            int subMatcher = indigoSessionPtr->_checkResult(bingoSearchSub(db, query.id(), ""));
            int resultObj = bingoGetObject(subMatcher);
            while (bingoNext(subMatcher) && result.size() < ansCount && !queryData.checkTimeOut() &&
                   !queryData.checkShouldStop()) {
                int clone = indigoClone(resultObj);
                result.emplace_back(clone);
            }
            if (queryData.checkTimeOut()) {
                LOG(INFO) << "Search stopped due to timeout";
            }
            bingoEndSearch(subMatcher);
            queryData.addAnswers(result);
            queryData.tagFinishTask();
        }
    }

    BingoNoSQLSearchData::BingoNoSQLSearchData(int db, IndigoSessionPtr indigoSessionPtr,
                                               TimeTicker &timeTicker, size_t ansCount, size_t threadsCount,
                                               double timeLimit) :
            SearchData(timeTicker, ansCount, threadsCount, timeLimit), db(db),
            indigoSessionPtr(std::move(indigoSessionPtr)) {}

    unique_ptr <QueryData<CIDType>>
    BingoNoSQLSearchData::search(const string &querySmiles, const PropertiesFilter::Bounds &) {
        LOG(INFO) << "Start search: " << querySmiles;
        unique_ptr<IndigoQueryMolecule> molecule = nullptr;
        try {
            molecule = make_unique<IndigoQueryMolecule>(indigoSessionPtr->loadQueryMolecule(querySmiles));
            molecule->aromatize();
        }
        catch (const IndigoException &e) {
            LOG(ERROR) << "Cannot parse smiles: " << querySmiles << " (" << e.what() << ")";
            return nullptr;
        }

        auto queryData = make_unique<QueryData<CIDType>>(ansCount, timeLimit, make_unique<AlwaysTrueFilter<CIDType>>());
        queryData->addTask(
                std::async(launch::async, searchInBingoDB, ref(*queryData), ansCount, cref(indigoSessionPtr), db,
                           std::move(*molecule)));
        return std::move(queryData);
    }

    BingoNoSQLSearchData::~BingoNoSQLSearchData() {
        bingoCloseDatabase(db);
    }
} // qtr