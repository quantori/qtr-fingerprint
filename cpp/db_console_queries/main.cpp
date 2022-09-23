#include <filesystem>
#include <string>

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
          "Path to file where splitter tree are stored");

ABSL_FLAG(std::string, db_name, "",
          "Name of folders with data base's files");

ABSL_FLAG(uint64_t, start_search_depth, -1,
          "Depth for start search from. There will be 2^start_search_depth threads.");

void initLogging(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    google::SetLogDestination(google::INFO, "db_console_queries.info");
    FLAGS_alsologtostderr = true;
}

struct Args {
    std::vector<std::filesystem::path> dataDirPaths;
    std::filesystem::path otherDataPath;
    std::string dbName;
    uint64_t startSearchDepth;

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

        for (auto &dir: dataDirPaths) {
            dbDataDirsPaths.emplace_back(dir / dbName);
        }
        for (size_t i = 0; i < dbDataDirsPaths.size(); i++) {
            LOG(INFO) << "dbDataDirPaths[" << i << "]: " << dbDataDirsPaths[i];
        }

        dbOtherDataPath = otherDataPath / dbName;
        LOG(INFO) << "dbOtherDataPath" << dbOtherDataPath;

        ballTreePath = dbOtherDataPath / "tree";
        LOG(INFO) << "splitterTreePath: " << ballTreePath;

        smilesTablePath = dbOtherDataPath / "smilesTable";
        LOG(INFO) << "smilesTablePath: " << smilesTablePath;

        fingerprintTablesPath = dbOtherDataPath / "fingerprintTables";
        LOG(INFO) << "fingerprintTablePaths" << fingerprintTablesPath;
    }
};

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

    size_t ansCount = 10;
//    std::cout << "Enter number of answers you want to find: ";
//    std::cin >> ansCount;

    std::string smiles;
    while (std::cin >> smiles) {
        timeTicker.tick();
        qtr::IndigoFingerprint fingerprint;
        try {
            fingerprint = qtr::IndigoFingerprintFromSmiles(smiles);
        }
        catch (std::exception &exception) {
            std::cout << "skip query:" << exception.what() << std::endl;
            continue;
        }
//        auto it = std::find(smilesTable.begin(), smilesTable.end(), smiles);
//        if (it != smilesTable.end()) {
//            uint64_t fpId = it - smilesTable.begin();
//            std::cout << "Check fingerprint " << fpId << std::endl;
//            for (const auto &filePath: qtr::findFiles(args.dbOtherDataPath / "fingerprintTables", ".ft")) {
//                for (const auto &[id, fp]: qtr::FingerprintTableReader(filePath)) {
//                    if (id == fpId) {
//                        std::cout << "same fingerprints?: " << (fp == fingerprint ? "YES" : "NO") << std::endl;
//                        if (!(fp == fingerprint)) {
//                            std::cout << "Diff bits: ";
//                            for (size_t i = 0; i < qtr::IndigoFingerprint::size(); i++) {
//                                if (fp[i] != fingerprint[i]) {
//                                    std::cout << i << ' ';
//                                }
//                            }
//                            std::cout << std::endl;
//                        }
//                    }
//                }
//            }
//            std::cout << "finish fingerprint check" << std::endl;
//        }
        auto ans = ballTree.search(fingerprint, -1, args.startSearchDepth);
        std::cout << "found answers: " << ans.size() << std::endl;
        size_t answersToPrint = std::min(ansCount, ans.size());
        for (size_t i = 0; i < answersToPrint; i++) {
            if (i == 0)
                std::cout << "[";
            std::cout << '\"' << smilesTable[ans[i]] << '\"';
            std::cout << (i + 1 == answersToPrint ? "]\n" : ", ");
        }
        timeTicker.tick("Search time");
    }

    return 0;
}