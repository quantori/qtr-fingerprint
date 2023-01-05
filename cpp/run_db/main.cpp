#include <string>

#include <glog/logging.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>

#include "Utils.h"
#include "BallTreeRAMSearchEngine.h"
#include "HuffmanCoder.h"
#include "SmilesTable.h"
#include "RunDbUtils.h"

using namespace std;
using namespace qtr;

#include "modes/web/WebMode.h"
#include "modes/InteractiveMode.h"
#include "modes/FromFileMode.h"


ABSL_FLAG(vector<string>, data_dir_paths, {},
          "Path to directories where data are stored");

ABSL_FLAG(string, other_data_path, "",
          "Path to directory with other files");

ABSL_FLAG(string, db_name, "",
          "Name of folders with data base's files");

ABSL_FLAG(uint64_t, start_search_depth, -1,
          "Depth for start search from. There will be 2^start_search_depth threads.");

ABSL_FLAG(string, mode, "",
          "possible modes: interactive, from_file, web");



ABSL_FLAG(string, input_file, "",
          "file to load test molecules from");

ABSL_FLAG(uint64_t, ans_count, -1,
          "how many answers program will find");

struct Args {

    enum class Mode {
        Interactive,
        FromFile,
        Web
    };

    vector<filesystem::path> dataDirPaths;
    filesystem::path otherDataPath;
    string dbName;
    uint64_t startSearchDepth;
    Mode mode;
    filesystem::path inputFile;
    uint64_t ansCount;

    vector<filesystem::path> dbDataDirsPaths;
    filesystem::path dbOtherDataPath;

    filesystem::path ballTreePath;
    filesystem::path smilesTablePath;
    filesystem::path huffmanCoderPath;
    filesystem::path idToStringDirPath;

    Args(int argc, char *argv[]) {
        absl::ParseCommandLine(argc, argv);

        vector<string> dataDirPathsStrings = absl::GetFlag(FLAGS_data_dir_paths);
        copy(dataDirPathsStrings.begin(), dataDirPathsStrings.end(), back_inserter(dataDirPaths));
        emptyArgument(dataDirPaths, "Please specify data_dir_paths option");
        for (size_t i = 0; i < dataDirPaths.size(); i++) {
            LOG(INFO) << "dataDirPaths[" << i << "]: " << dataDirPaths[i];
        }

        otherDataPath = absl::GetFlag(FLAGS_other_data_path);
        emptyArgument(otherDataPath, "Please specify other_data_path option");
        LOG(INFO) << "otherDataPath: " << otherDataPath;

        dbName = absl::GetFlag(FLAGS_db_name);
        emptyArgument(dbName, "Please specify db_name option");
        LOG(INFO) << "dbName: " << dbName;

        startSearchDepth = absl::GetFlag(FLAGS_start_search_depth);
        if (startSearchDepth == -1) {
            LOG(INFO) << "Please specify start_search_depth option";
            exit(-1);
        }
        LOG(INFO) << "startSearchDepth: " << startSearchDepth;

        string modeStr = absl::GetFlag(FLAGS_mode);
        emptyArgument(modeStr, "Please specify mode option");
        inputFile = absl::GetFlag(FLAGS_input_file);
        LOG(INFO) << "inputFile: " << inputFile;
        if (modeStr == "interactive") {
            mode = Mode::Interactive;
            LOG(INFO) << "mode: interactive";
        } else if (modeStr == "web") {
            mode = Mode::Web;
            LOG(INFO) << "mode: web server";
        } else if (modeStr == "from_file") {
            mode = Mode::FromFile;
            LOG(INFO) << "mode: fromFile";
            emptyArgument(inputFile, "Please specify input_file option");
        } else {
            LOG(ERROR) << "Bad mode option value";
            exit(-1);
        }

        ansCount = absl::GetFlag(FLAGS_ans_count);
        LOG(INFO) << "_stopAnswersNumber: " << ansCount;

        for (auto &dir: dataDirPaths) {
            dbDataDirsPaths.emplace_back(dir / dbName);
        }
        for (size_t i = 0; i < dbDataDirsPaths.size(); i++) {
            LOG(INFO) << "dbDataDirPaths[" << i << "]: " << dbDataDirsPaths[i];
        }

        dbOtherDataPath = otherDataPath / dbName;
        LOG(INFO) << "dbOtherDataPath" << dbOtherDataPath;

        ballTreePath = dbOtherDataPath / "tree";
        LOG(INFO) << "ballTreePath: " << ballTreePath;

        smilesTablePath = dbOtherDataPath / "smilesTable";
        LOG(INFO) << "smilesTablePath: " << smilesTablePath;

        huffmanCoderPath = dbOtherDataPath / "huffman";
        LOG(INFO) << "huffmanCoderPath: " << huffmanCoderPath;

        idToStringDirPath = dbOtherDataPath / "id_string";
        LOG(INFO) << "idToStringDirPath: " << idToStringDirPath;
    }
};

int main(int argc, char *argv[]) {
    initLogging(argv, google::INFO, "run_db.info", true);
    Args args(argc, argv);

    TimeTicker timeTicker;
    BufferedReader ballTreeReader(args.ballTreePath);

    HuffmanCoder huffmanCoder = HuffmanCoder::load(args.huffmanCoderPath);

    auto loadSmilesTableTask = async(launch::async, loadSmilesTable, cref(args.smilesTablePath), cref(huffmanCoder));
    LOG(INFO) << "Start ball tree loading";
    BallTreeRAMSearchEngine ballTree(ballTreeReader, args.dbDataDirsPaths);
    LOG(INFO) << "Finish ball tree loading";
    auto smilesTablePtr = loadSmilesTableTask.get();
    timeTicker.tick("DB initialization");
    RunMode *mode = nullptr;
    if (args.mode == Args::Mode::Interactive)
        mode = dynamic_cast<RunMode *>(new InteractiveMode(ballTree, smilesTablePtr, timeTicker, args.ansCount,
                                                           args.startSearchDepth));
    else if (args.mode == Args::Mode::FromFile)
        mode = dynamic_cast<RunMode *>(new FromFileMode(ballTree, smilesTablePtr, timeTicker, args.inputFile,
                                                        args.ansCount,
                                                        args.startSearchDepth));
    else if (args.mode == Args::Mode::Web)
        mode = dynamic_cast<RunMode *>(new WebMode(ballTree, smilesTablePtr, timeTicker, args.ansCount,
                                                   args.startSearchDepth, args.idToStringDirPath));

    mode->run();

    return 0;
}