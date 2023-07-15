#include "RunQtr.h"
#include "Args.h"
#include "TimeTicker.h"
#include "modes/RunMode.h"
#include "modes/InteractiveMode.h"
#include "modes/FromFileMode.h"
#include "modes/web/WebMode.h"
#include "search_data/RamSearchData.h"
#include "search_data/DriveSearchData.h"
#include "properties_table_io/PropertiesTableReader.h"

using namespace std;

namespace qtr {
    namespace {
        std::shared_ptr<SmilesTable>
        loadSmilesTable(const filesystem::path &smilesTablePath, const HuffmanCoder &huffmanCoder) {
            LOG(INFO) << "Start smiles table loading";
            HuffmanSmilesTable::Builder builder(huffmanCoder);
            for (const auto &pair: StringTableReader(smilesTablePath)) {
                builder += pair;
            }
            LOG(INFO) << "Finish smiles table loading";
            return builder.buildPtr();
        }

        shared_ptr<BallTreeSearchEngine> loadBallTree(const Args &args) {
            BufferedReader ballTreeReader(args.ballTreePath());
            LOG(INFO) << "Start ball tree loading";
            shared_ptr<BallTreeSearchEngine> res;
            if (args.dbType() == Args::DataBaseType::QtrRam)
                res = make_shared<BallTreeRAMSearchEngine>(ballTreeReader, args.dbDataDirs());
            else {
                res = make_shared<BallTreeDriveSearchEngine>(ballTreeReader, args.dbDataDirs());
            }
            LOG(INFO) << "Finish ball tree loading";
            return res;
        }

        shared_ptr<IdConverter> loadIdConverter(const std::filesystem::path &idToStringDirPath) {
            return std::make_shared<IdConverter>(idToStringDirPath);
        }

        shared_ptr<vector<PropertiesFilter::Properties>> loadPropertiesTable(const std::filesystem::path &propertiesTablePath) {
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

        shared_ptr<SearchData> loadRamSearchData(const Args &args, TimeTicker &timeTicker) {
            HuffmanCoder huffmanCoder = HuffmanCoder::load(args.huffmanCoderPath());
            auto loadBallTreeTask = async(launch::async, loadBallTree, cref(args));
            auto loadSmilesTableTask = async(launch::async, loadSmilesTable, args.smilesTablePath(), cref(huffmanCoder));
            auto loadIdConverterTask = async(launch::async, loadIdConverter, args.idToStringDir());
            auto loadPropertyTableTask = async(launch::async, loadPropertiesTable, args.propertyTablePath());

            auto ballTreePtr = loadBallTreeTask.get();
            auto smilesTablePtr = loadSmilesTableTask.get();
            auto idConverterPtr = loadIdConverterTask.get();
            auto propertiesTablePtr = loadPropertyTableTask.get();

            return make_shared<RamSearchData>(ballTreePtr, idConverterPtr, timeTicker, args.ansCount(), args.threads(),
                                              args.timeLimit(), smilesTablePtr, propertiesTablePtr);
        }

        shared_ptr<SearchData> loadDriveSearchData(const Args &args, TimeTicker &timeTicker) {
            auto loadBallTreeTask = async(launch::async, loadBallTree, cref(args));
            auto loadIdConverterTask = async(launch::async, loadIdConverter, args.idToStringDir());

            auto ballTreePtr = loadBallTreeTask.get();
            auto idConverterPtr = loadIdConverterTask.get();

            return make_shared<DriveSearchData>(ballTreePtr, idConverterPtr, timeTicker, args.ansCount(), args.threads(),
                                                args.timeLimit());
        }
        
        shared_ptr<SearchData> loadSearchData(const Args &args, TimeTicker &timeTicker) {
            if (args.dbType() == Args::DataBaseType::QtrRam) {
                return loadRamSearchData(args, timeTicker);
            } else if (args.dbType() == Args::DataBaseType::QtrDrive) {
                return loadDriveSearchData(args, timeTicker);
            }
            throw std::logic_error("Undefined db type");
        }
    }
    
    void runQtrDB(const Args& args) {
        TimeTicker timeTicker;
        auto searchData = loadSearchData(args, timeTicker);
        timeTicker.tick("Db data loading");

        shared_ptr<RunMode> mode = nullptr;
        if (args.mode() == Args::Mode::Interactive)
            mode = make_shared<InteractiveMode>(searchData);
        else if (args.mode() == Args::Mode::FromFile)
            mode = make_shared<FromFileMode>(searchData, args.queriesFile());
        else if (args.mode() == Args::Mode::Web)
            mode = make_shared<WebMode>(searchData);
        mode->run();
    }
}
