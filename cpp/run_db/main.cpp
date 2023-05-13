#include <string>
#include <stdexcept>

#include <glog/logging.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>

#include "Utils.h"
#include "BallTreeRAMSearchEngine.h"
#include "HuffmanCoder.h"
#include "SmilesTable.h"
#include "RunDbUtils.h"
#include "properties_table_io/PropertiesTableReader.h"
#include "search_data/RamSearchData.h"
#include "search_data/DriveSearchData.h"
#include "modes/web/WebMode.h"
#include "modes/InteractiveMode.h"
#include "modes/FromFileMode.h"


using namespace std;
using namespace qtr;


ABSL_FLAG(vector<string>, data_dir_paths, {},
          "Path to directories where data are stored");

ABSL_FLAG(string, other_data_path, "",
          "Path to directory with other files");

ABSL_FLAG(string, db_name, "",
          "Name of folders with data base's files");

ABSL_FLAG(uint64_t, threads_count, -1,
          "Number of threads to process leafs.");

ABSL_FLAG(string, mode, "",
          "possible modes: interactive, from_file, web");

ABSL_FLAG(string, input_file, "",
          "file to load test molecules from");

ABSL_FLAG(uint64_t, ans_count, -1,
          "how many answers program will find");

ABSL_FLAG(string, db_type, "",
          "Possible types: on_drive, in_ram");

struct Args {

    enum class Mode {
        Interactive,
        FromFile,
        Web
    };

    vector<filesystem::path> destDirPaths;
    filesystem::path otherDataPath;
    string dbName;
    uint64_t threadsCount;
    Mode mode;
    filesystem::path inputFile;
    uint64_t ansCount;
    DbType dbType;

    vector<filesystem::path> dbDataDirsPaths;
    filesystem::path dbOtherDataPath;

    filesystem::path ballTreePath;
    filesystem::path smilesTablePath;
    filesystem::path huffmanCoderPath;
    filesystem::path idToStringDirPath;
    filesystem::path propertyTableDestinationPath;

    Args(int argc, char *argv[]) {
        absl::ParseCommandLine(argc, argv);

        // get flags
        ansCount = absl::GetFlag(FLAGS_ans_count);
        vector<string> dataDirPathsStrings = absl::GetFlag(FLAGS_data_dir_paths);
        otherDataPath = absl::GetFlag(FLAGS_other_data_path);
        dbName = absl::GetFlag(FLAGS_db_name);
        threadsCount = absl::GetFlag(FLAGS_threads_count);
        string modeStr = absl::GetFlag(FLAGS_mode);
        inputFile = absl::GetFlag(FLAGS_input_file);
        string dbTypeStr = absl::GetFlag(FLAGS_db_type);

        // check empty flags
        checkEmptyArgument(inputFile, "Please specify input_file option");
        checkEmptyArgument(dataDirPathsStrings, "Please specify data_dir_paths option");
        checkEmptyArgument(otherDataPath, "Please specify other_data_path option");
        checkEmptyArgument(dbName, "Please specify db_name option");
        checkEmptyArgument(modeStr, "Please specify mode option");
        checkEmptyArgument(dbTypeStr, "Please specify db_type option");
        if (threadsCount == -1) {
            LOG(INFO) << "Please specify threads_count option";
            exit(-1);
        }

        // init values
        copy(dataDirPathsStrings.begin(), dataDirPathsStrings.end(), back_inserter(destDirPaths));
        for (auto &dir: destDirPaths) {
            dbDataDirsPaths.emplace_back(dir / dbName);
        }
        dbOtherDataPath = otherDataPath / dbName;
        ballTreePath = dbOtherDataPath / "tree";
        smilesTablePath = dbOtherDataPath / "smilesTable";
        huffmanCoderPath = dbOtherDataPath / "huffman";
        idToStringDirPath = dbOtherDataPath / "idToString";
        propertyTableDestinationPath = dbOtherDataPath / "propertyTable";
        if (modeStr == "interactive") {
            mode = Mode::Interactive;
        } else if (modeStr == "web") {
            mode = Mode::Web;
        } else if (modeStr == "from_file") {
            mode = Mode::FromFile;
        } else {
            LOG(ERROR) << "Bad mode option value";
            exit(-1);
        }
        if (dbTypeStr == "in_ram") {
            dbType = DbType::InRam;
        } else if (dbTypeStr == "on_drive") {
            dbType = DbType::OnDrive;
        } else {
            LOG(ERROR) << "Bad mode option value";
            exit(-1);
        }

        // log
        LOG(INFO) << "inputFile: " << inputFile;
        LOG(INFO) << "otherDataPath: " << otherDataPath;
        LOG(INFO) << "dbName: " << dbName;
        LOG(INFO) << "threadsCount: " << threadsCount;
        LOG(INFO) << "dbOtherDataPath" << dbOtherDataPath;
        LOG(INFO) << "ballTreePath: " << ballTreePath;
        LOG(INFO) << "smilesTablePath: " << smilesTablePath;
        LOG(INFO) << "huffmanCoderPath: " << huffmanCoderPath;
        LOG(INFO) << "idToStringDirPath: " << idToStringDirPath;
        LOG(INFO) << "propertyTableDestinationPath: " << propertyTableDestinationPath;
        for (size_t i = 0; i < dbDataDirsPaths.size(); i++) {
            LOG(INFO) << "dbDataDirPaths[" << i << "]: " << dbDataDirsPaths[i];
        }
        LOG(INFO) << "stopAnswersNumber: " << ansCount;
        for (size_t i = 0; i < destDirPaths.size(); i++) {
            LOG(INFO) << "destDirPaths[" << i << "]: " << destDirPaths[i];
        }
        if (mode == Mode::FromFile) {
            LOG(INFO) << "mode: from file";
        } else if (mode == Mode::Interactive) {
            LOG(INFO) << "mode: interactive";
        } else if (mode == Mode::Web) {
            LOG(INFO) << "mode: web server";
        }
    }
};

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
    BufferedReader ballTreeReader(args.ballTreePath);
    LOG(INFO) << "Start ball tree loading";
    shared_ptr<BallTreeSearchEngine> res;
    if (args.dbType == DbType::InRam)
        res = make_shared<BallTreeRAMSearchEngine>(ballTreeReader, args.dbDataDirsPaths);
    else {
        res = make_shared<BallTreeDriveSearchEngine>(ballTreeReader, args.dbDataDirsPaths);
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
    HuffmanCoder huffmanCoder = HuffmanCoder::load(args.huffmanCoderPath);
    auto loadBallTreeTask = async(launch::async, loadBallTree, cref(args));
    auto loadSmilesTableTask = async(launch::async, loadSmilesTable, cref(args.smilesTablePath), cref(huffmanCoder));
    auto loadIdConverterTask = async(launch::async, loadIdConverter, cref(args.idToStringDirPath));
    auto loadPropertyTableTask = async(launch::async, loadPropertiesTable, cref(args.propertyTableDestinationPath));

    auto ballTreePtr = loadBallTreeTask.get();
    auto smilesTablePtr = loadSmilesTableTask.get();
    auto idConverterPtr = loadIdConverterTask.get();
    auto propertiesTablePtr = loadPropertyTableTask.get();

    return make_shared<RamSearchData>(ballTreePtr, idConverterPtr, timeTicker, args.ansCount, args.threadsCount,
                                      smilesTablePtr, propertiesTablePtr);
}

shared_ptr<SearchData> loadDriveSearchData(const Args &args, TimeTicker &timeTicker) {
    auto loadBallTreeTask = async(launch::async, loadBallTree, cref(args));
    auto loadIdConverterTask = async(launch::async, loadIdConverter, cref(args.idToStringDirPath));

    auto ballTreePtr = loadBallTreeTask.get();
    auto idConverterPtr = loadIdConverterTask.get();

    return make_shared<DriveSearchData>(ballTreePtr, idConverterPtr, timeTicker, args.ansCount, args.threadsCount);
}

shared_ptr<SearchData> loadSearchData(const Args &args, TimeTicker &timeTicker) {
    if (args.dbType == qtr::DbType::InRam) {
        return loadRamSearchData(args, timeTicker);
    } else if (args.dbType == qtr::DbType::OnDrive) {
        return loadDriveSearchData(args, timeTicker);
    }
    throw std::logic_error("Undefined db type");
}

void runDb(const Args &args) {

    TimeTicker timeTicker;
    auto searchData = loadSearchData(args, timeTicker);
    timeTicker.tick("Db data loading");

    shared_ptr<RunMode> mode = nullptr;
    if (args.mode == Args::Mode::Interactive)
        mode = make_shared<InteractiveMode>(searchData);
    else if (args.mode == Args::Mode::FromFile)
        mode = make_shared<FromFileMode>(searchData, args.inputFile);
    else if (args.mode == Args::Mode::Web)
        mode = make_shared<WebMode>(searchData);
    mode->run();

}

int main(int argc, char *argv[]) {
    initLogging(argv, google::INFO, "run_db.info", true);
    Args args(argc, argv);

    runDb(args);

    return 0;
}