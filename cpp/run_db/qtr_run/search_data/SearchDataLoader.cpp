#include "SearchDataLoader.h"
#include "RunArgs.h"
#include "HuffmanCoder.h"

#include "QtrRamSearchData.h"
#include "QtrDriveSearchData.h"
#include "BingoNoSQLSearchData.h"
#include "RDKitSearchData.h"
#include "QtrEnumerationSearchData.h"
#include "properties_table_io/PropertiesTableReader.h"
#include "BallTreeRAMSearchEngine.h"
#include "HuffmanSmilesTable.h"
#include "string_table_io/StringTableReader.h"
#include "Profiling.h"

#include "IndigoException.h"
#include "base_cpp/scanner.h"
#include "src/bingo_object.h"
#include "molecule/smiles_loader.h"

#include "GraphMol/SubstructLibrary/SubstructLibrary.h"

using namespace std;
using namespace indigo_cpp;
using namespace indigo;
using namespace bingo;
using namespace qtr;

namespace {
    shared_ptr<SmilesTable>
    loadSmilesTable(const filesystem::path &smilesTablePath, const HuffmanCoder &huffmanCoder) {
        LOG(INFO) << "Start smiles table loading";
        HuffmanSmilesTable::Builder builder(huffmanCoder);
        for (const auto &pair: StringTableReader(smilesTablePath)) {
            builder += pair;
        }
        LOG(INFO) << "Finish smiles table loading";
        return builder.buildPtr();
    }

    shared_ptr<CFStorage> loadCFStorage(const filesystem::path &smilesTablePath) {
        LOG(INFO) << "Start CFStorage table loading";
        auto result = make_shared<CFStorage>();
        size_t counter = 0;
        for (const auto &[id, smiles]: StringTableReader(smilesTablePath)) {
            try {
                BufferScanner scanner(smiles.c_str(), smiles.size(), false);
                SmilesLoader loader(scanner);
                Molecule molecule;
                loader.loadMolecule(molecule);
                molecule.aromatize(AromaticityOptions());
                bingo::IndexMolecule indexMolecule(molecule, AromaticityOptions());
                auto &cfArr = result->Add(id, std::move(Array<char>()));
                indexMolecule.buildCfString(cfArr);
                counter++;
            }
            catch (const std::exception &e) {
                LOG_ERROR_AND_EXIT(e.what());
            }
            if (counter % 100000 == 0) {
                LOG(INFO) << "CFStorage loading: processed " << counter << " molecules";
            }
        }
        LOG(INFO) << "Finish CFStorage table loading";
        return result;
    }

    size_t loadFingerprintLength(const RunArgs &args) {
        ifstream in(args.fingerprintLengthFile());
        size_t res;
        in >> res;
        return res;
    }

    size_t loadMoleculesCount(const RunArgs &args) {
        ifstream in(args.totalMoleculesFile());
        size_t res;
        in >> res;
        return res;
    }

    shared_ptr<BallTreeSearchEngine> loadBallTree(const RunArgs &args) {
        BufferedReader ballTreeReader(args.ballTreePath());
        LOG(INFO) << "Start ball tree loading";
        shared_ptr<BallTreeSearchEngine> res;
        size_t fingerprintLength = loadFingerprintLength(args);
        size_t moleculesCount = loadMoleculesCount(args);
        if (args.dbType() == DatabaseType::QtrRam)
            res = make_shared<BallTreeRAMSearchEngine>(ballTreeReader, args.dbDataDirs(), fingerprintLength,
                                                       moleculesCount);
        else {
            res = make_shared<BallTreeDriveSearchEngine>(ballTreeReader, args.dbDataDirs(), fingerprintLength,
                                                         moleculesCount);
        }
        LOG(INFO) << "Finish ball tree loading";
        return res;
    }

    shared_ptr<vector<Fingerprint>>
    ballTreeToArray(const shared_ptr<BallTreeSearchEngine> &ballTree, size_t totalMolecules) {
        auto result = make_shared<vector<Fingerprint>>(totalMolecules);
        for (size_t leafId: ballTree->getLeafIds()) {
            for (auto &[i, fingerprint]: ballTree->getLeafContent(leafId)) {
                assert(i < result->size());
                result->at(i) = fingerprint;
            }
        }
        return result;
    }

    shared_ptr<IdConverter> loadIdConverter(const filesystem::path &idToStringDirPath) {
        return make_shared<IdConverter>(idToStringDirPath);
    }

    shared_ptr<vector<PropertiesFilter::Properties>>
    loadPropertiesTable(const filesystem::path &propertiesTablePath) {
        LOG(INFO) << "Start properties table loading";
        auto res = make_shared<vector<PropertiesFilter::Properties>>();
        auto reader = PropertiesTableReader(propertiesTablePath);
        for (const auto &[id, properties]: reader) {
            assert(id == res->size());
            res->emplace_back(properties);
        }
        LOG(INFO) << "Finish properties table loading";
        return res;
    }

    shared_ptr<SearchData> loadQtrRamSearchData(const RunArgs &args) {
        HuffmanCoder huffmanCoder = HuffmanCoder::load(args.huffmanCoderPath());
        auto loadBallTreeTask = async(launch::async, loadBallTree, cref(args));
        auto loadCFStorageTask = [&args]() {
            if (args.verificationStage())
                return async(launch::async, loadCFStorage, args.smilesTablePath());
            else
                return async(launch::async, []() -> shared_ptr<CFStorage> { return nullptr; });
        }();
        auto loadIdConverterTask = async(launch::async, loadIdConverter, args.idToStringDir());

        shared_ptr<vector<PropertiesFilter::Properties>> propertiesTablePtr = nullptr;
        if (args.properties()) {
            auto loadPropertyTableTask = async(launch::async, loadPropertiesTable, args.propertyTablePath());
            propertiesTablePtr = loadPropertyTableTask.get();
        }
        auto ballTreePtr = loadBallTreeTask.get();
        auto cfStoragePtr = loadCFStorageTask.get();
        auto idConverterPtr = loadIdConverterTask.get();

        return make_shared<QtrRamSearchData>(ballTreePtr, idConverterPtr, args.ansCount(), args.threads(),
                                             args.timeLimit(), cfStoragePtr, propertiesTablePtr,
                                             args.verificationStage());
    }

    shared_ptr<SearchData> loadQtrEnumerationSearchData(const RunArgs &args) {
        auto searchData = dynamic_pointer_cast<QtrRamSearchData>(loadQtrRamSearchData(args));
        assert(searchData != nullptr);
        auto ans = make_shared<QtrEnumerationSearchData>(*searchData);
        return ans;
    }

    shared_ptr<SearchData> loadQtrDriveSearchData(const RunArgs &args) {
        auto loadBallTreeTask = async(launch::async, loadBallTree, cref(args));
        auto loadIdConverterTask = async(launch::async, loadIdConverter, args.idToStringDir());

        auto ballTreePtr = loadBallTreeTask.get();
        auto idConverterPtr = loadIdConverterTask.get();

        return make_shared<QtrDriveSearchData>(ballTreePtr, idConverterPtr, args.ansCount(), args.threads(),
                                               args.timeLimit(), args.verificationStage());
    }

    shared_ptr<SearchData> loadBingoNoSQLSearchData(const RunArgs &args) {
        filesystem::path dbDataDir = args.dbDataDirs()[0];
        try {
            return make_shared<BingoNoSQLSearchData>(dbDataDir, args.ansCount(), args.threads(), args.timeLimit(),
                                                     args.verificationStage());
        }
        catch (const IndigoException &e) {
            LOG_ERROR_AND_EXIT(e.what());
        }
    }

    shared_ptr<SearchData> loadRDKitSearchData(const RunArgs &args) {
        filesystem::path dbDataDir = args.dbDataDirs()[0];
        filesystem::path moleculesDir = dbDataDir / "molecules";
        filesystem::path fingerprintsDir = dbDataDir / "fingerprintTables";
        try {
            return make_shared<RDKitSearchData>(moleculesDir, fingerprintsDir, args.ansCount(), args.threads(),
                                                args.timeLimit(), args.verificationStage());
        }
        catch (const exception &e) {
            LOG_ERROR_AND_EXIT(e.what());
        }
    }

}

shared_ptr<SearchData> SearchDataLoader::load(const RunArgs &args) {
    ProfileScope("SearchData loading");
    if (args.dbType() == DatabaseType::QtrRam) {
        return loadQtrRamSearchData(args);
    } else if (args.dbType() == DatabaseType::QtrDrive) {
        return loadQtrDriveSearchData(args);
    } else if (args.dbType() == DatabaseType::BingoNoSQL) {
        return loadBingoNoSQLSearchData(args);
    } else if (args.dbType() == DatabaseType::QtrEnumeration) {
        return loadQtrEnumerationSearchData(args);
    } else if (args.dbType() == DatabaseType::RDKit) {
        return loadRDKitSearchData(args);
    }
    throw logic_error("Undefined db type");
}
