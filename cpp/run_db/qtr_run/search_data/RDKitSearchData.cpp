#include "RDKitSearchData.h"
#include "data_io/fingerprint_table_io/FingerprintTableReader.h"
#include <GraphMol/Fingerprints/Fingerprints.h>

using namespace std;

namespace qtr {

    namespace {
        vector <pair<filesystem::path, filesystem::path>>
        zipDirs(const filesystem::path &dir1, const filesystem::path &dir2) {
            map<std::string, filesystem::path> dir2Files;
            for (auto &entry: filesystem::directory_iterator(dir2)) {
                if (!entry.is_regular_file()) {
                    continue;
                }
                dir2Files[entry.path().stem().string()] = entry.path();
            }
            vector<pair<filesystem::path, filesystem::path>> result;
            for (auto &entry: filesystem::directory_iterator(dir1)) {
                auto stem = entry.path().stem().string();
                auto it = dir2Files.find(stem);
                if (it == dir2Files.end()) {
                    throw runtime_error("Cannot zip directories " + dir1.string() + " and " + dir2.string() +
                                        " because their contents are not aligned. " + entry.path().filename().string() +
                                        " is presented only in the first dir");
                }
                result.emplace_back(entry.path(), it->second);
                dir2Files.erase(it);
            }
            if (!dir2Files.empty()) {
                throw runtime_error("Cannot zip directories " + dir1.string() + " and " + dir2.string() +
                                    " because their contents are not aligned. " +
                                    dir2Files.begin()->second.filename().string() +
                                    " is presented only in the second dir");
            }
            return result;
        }


        void
        searchInRDKitDB(QueryData<CIDType> &queryData, uint64_t ansCount, const RDKit::SubstructLibrary &substructLib,
                        const std::shared_ptr<RDKit::ROMol> &mol) {
            RDKit::SubstructMatchParameters params;
            params.recursionPossible = true;
            params.useChirality = true;
            params.useQueryQueryMatches = false;
            int maxResults = ansCount == uint64_t(-1) ? -1 : int(ansCount);
            const unsigned int STEP = 100; // check timeout every STEP iterations
            for (unsigned int block = 0; block < substructLib.size() && maxResults != 0; block += STEP) {
                auto matches = substructLib.getMatches(*mol, block, std::min(substructLib.size(), block + STEP), params,
                                                       1, maxResults);
                if (!matches.empty()) {
                    vector<CIDType> results;
                    results.reserve(matches.size());
                    for (auto &i: matches) {
                        results.emplace_back(stoi(substructLib.getKeys().getKey(i)));
                    }
                    queryData.addAnswers(results);
                }if (queryData.checkShouldStop() || queryData.checkTimeOut()) {
                    break;
                }
                if (maxResults != -1) {
                    maxResults -= (int) matches.size();
                }
            }
            if (queryData.checkTimeOut()) {
                LOG(INFO) << "Search stopped due to timeout";
            }
            queryData.tagFinishTask();
        }


        void addElementsToHandlers(const filesystem::path &moleculesFile, const filesystem::path &fingerprintTableFile,
                                   const boost::shared_ptr<RDKit::CachedMolHolder> &molHandler,
                                   const boost::shared_ptr<RDKit::PatternHolder> &patternHolder,
                                   const boost::shared_ptr<RDKit::KeyFromPropHolder> &keyHolder,
                                   std::mutex &mutex) {
            ifstream moleculesIn(moleculesFile, ios::binary);
            for (auto [id, fp]: FingerprintTableReader(fingerprintTableFile)) {
                RDKit::ROMol mol;
                try {
                    RDKit::MolPickler::molFromPickle(moleculesIn, mol);
                } catch (const exception &e) {
                    LOG_ERROR_AND_EXIT("Cannot extract molecule from pickle: "s + e.what());
                }
                auto bitVector = new ExplicitBitVect(fp.size());
                for (size_t i = 0; i < fp.size(); i++) {
                    if (fp[i]) {
                        bitVector->setBit(i);
                    }
                }
                {
                    lock_guard<std::mutex> guard(mutex);
                    molHandler->addMol(mol);
                    patternHolder->addFingerprint(bitVector);
                    keyHolder->addKey(to_string(id));
                }
            }

        }
    }


    RDKitSearchData::RDKitSearchData(const std::filesystem::path &moleculesDir,
                                     const std::filesystem::path &fingerprintTablesDir,
                                     size_t ansCount, size_t threadsCount, double timeLimit, bool verificationStage) :
            SearchData(ansCount, threadsCount, timeLimit, verificationStage) {
        auto molHandler = boost::make_shared<RDKit::CachedMolHolder>();
        auto fpHandler = boost::make_shared<RDKit::PatternHolder>();
        auto idHandler = boost::make_shared<RDKit::KeyFromPropHolder>();
        std::mutex mutex;
        vector<future<void>> tasks;
        for (const auto &[moleculesFile, fingerprintsFile]: zipDirs(moleculesDir, fingerprintTablesDir)) {
            tasks.push_back(async(launch::async,
                                  [moleculesFile, fingerprintsFile, &molHandler, &fpHandler, &idHandler, &mutex] {
                                      addElementsToHandlers(moleculesFile, fingerprintsFile, molHandler, fpHandler,
                                                            idHandler, mutex);
                                  }));
        }
        for (auto &task: tasks) {
            task.wait();
        }
        _substructLibrary = std::make_shared<RDKit::SubstructLibrary>(molHandler, fpHandler, idHandler);
    }

    unique_ptr <QueryData<CIDType>>
    RDKitSearchData::search(const SearchData::Query &query, const PropertiesFilter::Bounds &) {
        assert(query.smiles != nullptr);
        assert(query.fingerprint == nullptr);
        LOG(INFO) << "Start search: " << *query.smiles;
        shared_ptr<RDKit::ROMol> mol = nullptr;
        try {
            mol = std::shared_ptr<RDKit::ROMol>(RDKit::SmilesToMol(*query.smiles));
            if (mol == nullptr) {
                return nullptr;
            }
        } catch (const exception &e) {
            LOG(ERROR) << e.what();
            return nullptr;
        }

        auto queryData = make_unique<QueryData<
                CIDType >>(ansCount, timeLimit, make_unique<AlwaysTrueFilter<CIDType >>());
        queryData->addTask(
                async(launch::async, searchInRDKitDB, ref(*queryData), ansCount, cref(*_substructLibrary), mol));
        return queryData;
    }
} // qtr