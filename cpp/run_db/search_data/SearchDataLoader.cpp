#include "SearchDataLoader.h"
#include "Args.h"
#include "HuffmanCoder.h"

#include "QtrRamSearchData.h"
#include "QtrDriveSearchData.h"
#include "BingoNoSQLSearchData.h"
#include "properties_table_io/PropertiesTableReader.h"
#include "BallTreeRAMSearchEngine.h"
#include "HuffmanSmilesTable.h"
#include "string_table_io/StringTableReader.h"

#include "IndigoException.h"
#include "base_cpp/scanner.h"
#include "src/bingo_object.h"
#include "molecule/smiles_loader.h"

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
                logErrorAndExit(e.what());
            }
            if (counter % 100000 == 0) {
                LOG(INFO) << "CFStorage loading: processed " << counter << " molecules";
            }
        }
        LOG(INFO) << "Finish CFStorage table loading";
        return result;
    }

    shared_ptr<BallTreeSearchEngine> loadBallTree(const Args &args) {
        BufferedReader ballTreeReader(args.ballTreePath());
        LOG(INFO) << "Start ball tree loading";
        shared_ptr<BallTreeSearchEngine> res;
        if (args.dbType() == DatabaseType::QtrRam)
            res = make_shared<BallTreeRAMSearchEngine>(ballTreeReader, args.dbDataDirs());
        else {
            res = make_shared<BallTreeDriveSearchEngine>(ballTreeReader, args.dbDataDirs());
        }
        LOG(INFO) << "Finish ball tree loading";
        return res;
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

    shared_ptr<SearchData> loadQtrRamSearchData(const Args &args, TimeTicker &timeTicker) {
        HuffmanCoder huffmanCoder = HuffmanCoder::load(args.huffmanCoderPath());
        auto loadBallTreeTask = async(launch::async, loadBallTree, cref(args));
//        auto loadSmilesTableTask = async(launch::async, loadSmilesTable, args.smilesTablePath(),
//                                         cref(huffmanCoder));
        auto loadCFStorageTask = async(launch::async, loadCFStorage, args.smilesTablePath());
        auto loadIdConverterTask = async(launch::async, loadIdConverter, args.idToStringDir());
        auto loadPropertyTableTask = async(launch::async, loadPropertiesTable, args.propertyTablePath());

        auto ballTreePtr = loadBallTreeTask.get();
//        auto smilesTablePtr = loadSmilesTableTask.get();
        auto cfStoragePtr = loadCFStorageTask.get();
        auto idConverterPtr = loadIdConverterTask.get();
        auto propertiesTablePtr = loadPropertyTableTask.get();

        return make_shared<QtrRamSearchData>(ballTreePtr, idConverterPtr, timeTicker, args.ansCount(),
                                             args.threads(), args.timeLimit(), cfStoragePtr, propertiesTablePtr);
    }

    shared_ptr<SearchData> loadQtrDriveSearchData(const Args &args, TimeTicker &timeTicker) {
        auto loadBallTreeTask = async(launch::async, loadBallTree, cref(args));
        auto loadIdConverterTask = async(launch::async, loadIdConverter, args.idToStringDir());

        auto ballTreePtr = loadBallTreeTask.get();
        auto idConverterPtr = loadIdConverterTask.get();

        return make_shared<QtrDriveSearchData>(ballTreePtr, idConverterPtr, timeTicker, args.ansCount(),
                                               args.threads(), args.timeLimit());
    }

    shared_ptr<SearchData> loadBingoNoSQLSearchData(const Args &args, TimeTicker &timeTicker) {
        filesystem::path dbDataDir = args.dbDataDirs()[0];
        try {
            return make_shared<BingoNoSQLSearchData>(dbDataDir, timeTicker, args.ansCount(), args.threads(),
                                                     args.timeLimit());
        }
        catch (const IndigoException &e) {
            logErrorAndExit(e.what());
        }

    }
}

shared_ptr<SearchData> SearchDataLoader::load(const Args &args, TimeTicker &timeTicker) {
    if (args.dbType() == DatabaseType::QtrRam) {
        return loadQtrRamSearchData(args, timeTicker);
    } else if (args.dbType() == DatabaseType::QtrDrive) {
        return loadQtrDriveSearchData(args, timeTicker);
    } else if (args.dbType() == DatabaseType::BingoNoSQL) {
        return loadBingoNoSQLSearchData(args, timeTicker);
    }
    throw logic_error("Undefined db type");
}