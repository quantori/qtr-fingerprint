#include <filesystem>
#include <string>
#include <algorithm>
#include <numeric>

#include <glog/logging.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>

#include "Utils.h"
#include "io/BufferedReader.h"
#include "data_io/smiles_table_io/SmilesTableReader.h"
#include "BallTree.h"
#include "fingerprint_table_io/FingerprintTableReader.h"

ABSL_FLAG(std::vector<std::string>, data_dir_paths, {},
          "Path to directories where data are stored");

ABSL_FLAG(std::string, other_data_path, "",
          "Path to directory with other files");

ABSL_FLAG(std::string, db_name, "",
          "Name of folders with data base's files");

ABSL_FLAG(uint64_t, start_search_depth, -1,
          "Depth for start search from. There will be 2^start_search_depth threads.");

ABSL_FLAG(std::string, mode, "",
          "possible modes: interactive, from_file");

ABSL_FLAG(std::string, input_file, "",
          "file to load test molecules from");

ABSL_FLAG(std::uint64_t, ans_count, -1,
          "how many answers program will find");

void initLogging(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    google::SetLogDestination(google::INFO, "run_db.info");
    FLAGS_alsologtostderr = true;
}

struct Args {

    enum class Mode {
        Interactive,
        FromFile
    };

    std::vector<std::filesystem::path> dataDirPaths;
    std::filesystem::path otherDataPath;
    std::string dbName;
    uint64_t startSearchDepth;
    Mode mode;
    std::filesystem::path inputFile;
    uint64_t ansCount;

    std::vector<std::filesystem::path> dbDataDirsPaths;
    std::filesystem::path dbOtherDataPath;

    std::filesystem::path ballTreePath;
    std::filesystem::path smilesTablePath;
    std::filesystem::path fingerprintTablesPath;

    Args(int argc, char *argv[]) {
        absl::ParseCommandLine(argc, argv);

        std::vector<std::string> dataDirPathsStrings = absl::GetFlag(FLAGS_data_dir_paths);
        std::copy(dataDirPathsStrings.begin(), dataDirPathsStrings.end(), std::back_inserter(dataDirPaths));
        qtr::emptyArgument(dataDirPaths, "Please specify data_dir_paths option");
        for (size_t i = 0; i < dataDirPaths.size(); i++) {
            LOG(INFO) << "dataDirPaths[" << i << "]: " << dataDirPaths[i];
        }

        otherDataPath = absl::GetFlag(FLAGS_other_data_path);
        qtr::emptyArgument(otherDataPath, "Please specify other_data_path option");
        LOG(INFO) << "otherDataPath: " << otherDataPath;

        dbName = absl::GetFlag(FLAGS_db_name);
        qtr::emptyArgument(dbName, "Please specify db_name option");
        LOG(INFO) << "dbName: " << dbName;

        startSearchDepth = absl::GetFlag(FLAGS_start_search_depth);
        if (startSearchDepth == -1) {
            LOG(INFO) << "Please specify start_search_depth option";
            exit(-1);
        }
        LOG(INFO) << "startSearchDepth: " << startSearchDepth;

        std::string modeStr = absl::GetFlag(FLAGS_mode);
        qtr::emptyArgument(modeStr, "Please specify mode option");
        inputFile = absl::GetFlag(FLAGS_input_file);
        LOG(INFO) << "inputFile: " << inputFile;
        if (modeStr == "interactive") {
            mode = Mode::Interactive;
            LOG(INFO) << "mode: interactive";
        } else if (modeStr == "from_file") {
            mode = Mode::FromFile;
            LOG(INFO) << "mode: fromFile";
            qtr::emptyArgument(inputFile, "Please specify input_file option");
        } else {
            LOG(ERROR) << "Bad mode option value";
            exit(-1);
        }

        ansCount = absl::GetFlag(FLAGS_ans_count);
        LOG(INFO) << "ansCount: " << ansCount;

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

        fingerprintTablesPath = dbOtherDataPath / "fingerprintTables";
        LOG(INFO) << "fingerprintTablePaths" << fingerprintTablesPath;
    }
};

void runInteractive(const qtr::BallTree &ballTree, const std::vector<std::string> &smilesTable,
                    qtr::TimeTicker &timeTicker, const Args &args) {
    while (true) {
        std::cout << "Enter smiles: ";
        std::string smiles;
        std::cin >> smiles;
        if (smiles.empty())
            break;
        timeTicker.tick();
        qtr::IndigoFingerprint fingerprint;
        try {
            fingerprint = qtr::IndigoFingerprintFromSmiles(smiles);
        }
        catch (std::exception &exception) {
            std::cout << "skip query:" << exception.what() << std::endl;
            continue;
        }

        auto ans = ballTree.search(fingerprint, args.ansCount, args.startSearchDepth);
        std::cout << "found answers: " << ans.size() << std::endl;
        size_t answersToPrint = std::min(size_t(10), ans.size());
        for (size_t i = 0; i < answersToPrint; i++) {
            if (i == 0)
                std::cout << "[";
            std::cout << '\"' << smilesTable[ans[i]] << '\"';
            std::cout << (i + 1 == answersToPrint ? "]\n" : ", ");
        }
        std::cout<< "Search time: " << timeTicker.tick("Search time") << std::endl;
    }
}

void runFromFile(const qtr::BallTree &ballTree, const std::vector<std::string> &smilesTable,
                 qtr::TimeTicker &timeTicker, const Args &args) {
    std::ifstream input(args.inputFile);
    std::vector<std::string> queries;
    while (input.peek() != EOF) {
        std::string query;
        input >> query;
        std::string otherInfoInLine;
        std::getline(input, otherInfoInLine);
        queries.emplace_back(query);
    }
    std::vector<double> times;
    size_t skipped;
    LOG(INFO) << "Loaded " << queries.size() << " queries";
    for (size_t i = 0; i < queries.size(); i++) {
        LOG(INFO) << "Start search for " << i << ": " << queries[i];
        timeTicker.tick();
        qtr::IndigoFingerprint fingerprint;
        try {
            fingerprint = qtr::IndigoFingerprintFromSmiles(queries[i]);
        }
        catch (std::exception &exception) {
            std::cout << "skip query:" << exception.what() << std::endl;
            skipped++;
            continue;
        }
        auto ans = ballTree.search(fingerprint, args.ansCount, args.startSearchDepth);
        times.emplace_back(timeTicker.tick("search molecule " + std::to_string(i) + ": " + queries[i]));
    }

    double mean = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
    double min = *std::min_element(times.begin(), times.end());
    double max = *std::max_element(times.begin(), times.end());
    std::sort(times.begin(), times.end());
    double median = times.size() & 1 ? times[times.size() / 2] : times[times.size() / 2] + times[times.size() / 2 + 1];

    LOG(INFO) << "  mean: " << mean;
    LOG(INFO) << "   max: " << max;
    LOG(INFO) << "   min: " << min;
    LOG(INFO) << "median: " << median;

    std::cout << "  mean: " << mean << std::endl;
    std::cout << "   max: " << max << std::endl;
    std::cout << "   min: " << min << std::endl;
    std::cout << "median: " << median << std::endl;
}

int main(int argc, char *argv[]) {
    initLogging(argc, argv);
    Args args(argc, argv);

    qtr::TimeTicker timeTicker;
    std::vector<std::string> smilesTable;
    for (const auto &[id, smiles]: qtr::SmilesTableReader(args.smilesTablePath)) {
        if (id != smilesTable.size())
            throw std::exception();
        smilesTable.emplace_back(smiles);
    }
    qtr::BufferedReader ballTreeReader(args.ballTreePath);
    qtr::BallTree ballTree(ballTreeReader, args.dbDataDirsPaths);
    timeTicker.tick("DB initialization");

    if (args.mode == Args::Mode::Interactive)
        runInteractive(ballTree, smilesTable, timeTicker, args);
    else if (args.mode == Args::Mode::FromFile)
        runFromFile(ballTree, smilesTable, timeTicker, args);

    return 0;
}