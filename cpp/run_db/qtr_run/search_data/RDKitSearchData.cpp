#include "RDKitSearchData.h"

using namespace std;

namespace qtr {

    namespace {
        void searchInRDKitDB(QueryData<CIDType> &queryData, uint64_t ansCount, BingoMolecule &db,
                             const IndigoQueryMolecule &query) {
            // TODO: implement efficiently
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

    RDKitSearchData::RDKitSearchData(const std::filesystem::path &dbDataDir, size_t ansCount, size_t threadsCount,
                                     double timeLimit, bool verificationStage) :
            SearchData(ansCount, threadsCount, timeLimit, verificationStage) {
        // TODO

    }

    RDKitSearchData::~RDKitSearchData() {
        // TODO
    }

    unique_ptr <QueryData<CIDType>>
    RDKitSearchData::search(const SearchData::Query &query, const PropertiesFilter::Bounds &) {
        // TODO: rewrite for RDKit
        assert(query.smiles != nullptr);
        assert(query.fingerprint == nullptr);
        LOG(INFO) << "Start search: " << *query.smiles;
//        shared_ptr<RDKit::ROMol> molecule = nullptr;
        try {
//            molecule = shared_ptr<RDKit::ROMol>(RDKit::SmilesToMol(*query.smiles));
        }
        catch (const std::exception &e) {
//            LOG(ERROR) << "Cannot parse smiles: " << *query.smiles << " (" << e.what() << ")";
//            return nullptr;
        }
//
        auto queryData = make_unique<QueryData<CIDType>>(ansCount, timeLimit, make_unique<AlwaysTrueFilter<CIDType>>());
//        queryData->addTask(std::async(launch::async, searchInRDKitDB, ref(*queryData), ansCount, std::ref(db),
//                                      std::move(molecule)));
        return std::move(queryData);
    }
} // qtr