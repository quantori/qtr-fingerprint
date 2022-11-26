#include <filesystem>
#include <string>
#include <random>
#include <future>

#include <glog/logging.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>

#include "Utils.h"
#include "smiles_table_io/SmilesTableWriter.h"
#include "smiles_table_io/SmilesTableReader.h"
#include "fingerprint_table_io/FingerprintTableWriter.h"
#include "BallTreeBuilder.h"
#include "ball_tree/split_bit_selection/MaxDispersionBitSelector.h"

ABSL_FLAG(std::string, smiles_dir_path, "",
          "Path to directory where smiles tables are stored");

ABSL_FLAG(std::string, fingerprints_dir_path, "",
          "Path to directory where fingerprint tables are stored");

ABSL_FLAG(std::vector<std::string>, data_dir_paths, {},
          "Path to directories where data should be stored");

ABSL_FLAG(std::string, other_data_path, "",
          "Path to file where splitter tree should be stored");

ABSL_FLAG(uint64_t, subtree_parallel_depth, 0,
          "Depth on which subtree parallelization starts");

ABSL_FLAG(uint64_t, tree_depth, 0,
          "ball tree depth");

ABSL_FLAG(std::string, db_name, "",
          "Name of folders with data base's files");

namespace {
    std::mt19937 randomGenerator(0);
}

struct Args {
    std::filesystem::path smilesDirPath;
    std::filesystem::path fingerprintsDirPath;
    std::vector<std::filesystem::path> dataDirPaths;
    std::filesystem::path otherDataPath;
    uint64_t subtreeParallelDepth;
    uint64_t treeDepth;
    std::string dbName;

    std::vector<std::filesystem::path> dbDataDirsPaths;
    std::filesystem::path dbOtherDataPath;

    std::filesystem::path ballTreePath;
    std::filesystem::path smilesTablePath;

    Args(int argc, char *argv[]) {
        absl::ParseCommandLine(argc, argv);

        smilesDirPath = absl::GetFlag(FLAGS_smiles_dir_path);
        qtr::emptyArgument(smilesDirPath, "Please specify smiles_dir_path option");
        LOG(INFO) << "smilesDirPath: " << smilesDirPath;

        fingerprintsDirPath = absl::GetFlag(FLAGS_fingerprints_dir_path);
        qtr::emptyArgument(fingerprintsDirPath, "Please specify fingerprints_dir_path option");
        LOG(INFO) << "fingerprintsDirPath" << fingerprintsDirPath;

        std::vector<std::string> dataDirPathsStrings = absl::GetFlag(FLAGS_data_dir_paths);
        std::copy(dataDirPathsStrings.begin(), dataDirPathsStrings.end(), std::back_inserter(dataDirPaths));
        qtr::emptyArgument(dataDirPaths, "Please specify data_dir_paths option");
        for (size_t i = 0; i < dataDirPaths.size(); i++) {
            LOG(INFO) << "dataDirPaths[" << i << "]: " << dataDirPaths[i];
        }

        otherDataPath = absl::GetFlag(FLAGS_other_data_path);
        qtr::emptyArgument(otherDataPath, "Please specify other_data_path option");
        LOG(INFO) << "otherDataPath: " << otherDataPath;

        subtreeParallelDepth = absl::GetFlag(FLAGS_subtree_parallel_depth);
        qtr::emptyArgument(subtreeParallelDepth, "Please specify subtree_parallel_depth option");
        LOG(INFO) << "subtreeParallelDepth: " << subtreeParallelDepth;

        treeDepth = absl::GetFlag(FLAGS_tree_depth);
        qtr::emptyArgument(treeDepth, "Please specify tree_depth option");
        LOG(INFO) << "treeDepth: " << treeDepth;

        dbName = absl::GetFlag(FLAGS_db_name);
        if (dbName.empty()) {
            LOG(INFO) << "db_name option is not specified";
            dbName = qtr::generateDbName(dataDirPaths, otherDataPath);
            LOG(INFO) << "Generated name for data base: " << dbName;
        }
        LOG(INFO) << "dbName: " << dbName;

        for (auto &dir: dataDirPaths) {
            dbDataDirsPaths.emplace_back(dir / dbName);
        }
        for (size_t i = 0; i < dbDataDirsPaths.size(); i++) {
            LOG(INFO) << "dbDataDirPaths[" << i << "]: " << dbDataDirsPaths[i];
        }

        dbOtherDataPath = otherDataPath / dbName;
        LOG(INFO) << "dbOtherDataPath: " << dbOtherDataPath;

        ballTreePath = dbOtherDataPath / "tree";
        LOG(INFO) << "splitterTreePath: " << ballTreePath;

        smilesTablePath = dbOtherDataPath / "smilesTable";
        LOG(INFO) << "smilesTablePath: " << smilesTablePath;
    }
};

void initFileSystem(const Args &args) {
    std::vector<std::filesystem::path> alreadyExists;
    for (auto &dbDirPath: args.dbDataDirsPaths) {
        if (std::filesystem::exists(dbDirPath))
            alreadyExists.emplace_back(dbDirPath);
    }
    if (std::filesystem::exists(args.dbOtherDataPath)) {
        alreadyExists.emplace_back(args.dbOtherDataPath);
    }
    if (!alreadyExists.empty()) {
        std::cout << "Some data directories already exist: \n";
        for (auto &dir: alreadyExists) {
            std::cout << dir << '\n';
        }
        qtr::askAboutContinue("Data will be overridden");
        for (auto &dir: alreadyExists) {
            std::filesystem::remove_all(dir);
        }
    }
    for (auto &dbDirPath: args.dbDataDirsPaths) {
        std::filesystem::create_directory(dbDirPath);
        std::filesystem::create_directory(dbDirPath / "0");
    }
    std::filesystem::create_directory(args.dbOtherDataPath);
}

void distributeFingerprintTables(const Args &args) {
    std::vector<std::filesystem::path> ftFilePaths = qtr::findFiles(args.fingerprintsDirPath,
                                                                    qtr::fingerprintTableExtension);
    std::shuffle(ftFilePaths.begin(), ftFilePaths.end(), randomGenerator);
    size_t drivesCount = args.dbDataDirsPaths.size();
    for (size_t i = 0; i < ftFilePaths.size(); i++) {
        auto sourcePath = ftFilePaths[i];
        auto destinationPath = args.dbDataDirsPaths[i % drivesCount] / "0" / sourcePath.filename();
        std::filesystem::create_directory(destinationPath.parent_path());
        LOG(INFO) << "Copy " << sourcePath << " to " << destinationPath;
        std::filesystem::copy_file(sourcePath, destinationPath);
    }
}

size_t mergeSmilesTables(const Args &args) {
    std::vector<std::filesystem::path> stFilesPaths = qtr::findFiles(args.smilesDirPath, qtr::smilesTableExtension);
    std::vector<qtr::smiles_table_value_t> smilesTable;

    for (auto &stFile: stFilesPaths) {
        qtr::SmilesTableReader reader(stFile);
        std::copy(reader.begin(), reader.end(), std::back_inserter(smilesTable));
    }

    std::sort(smilesTable.begin(), smilesTable.end(),
              [](const qtr::smiles_table_value_t &a, const qtr::smiles_table_value_t &b) {
                  return a.first < b.first;
              });
    qtr::SmilesTableWriter writer(args.smilesTablePath);
    std::copy(smilesTable.begin(), smilesTable.end(), writer.begin());
    return smilesTable.size();
}

int main(int argc, char *argv[]) {
    qtr::initLogging(argv, google::INFO, "build_db.log", true);
    Args args(argc, argv);

    qtr::TimeTicker timeTicker;

    initFileSystem(args);
    timeTicker.tick("Filesystem initialization");
    auto smilesTask = std::async(std::launch::async, mergeSmilesTables, std::cref(args));
    auto fingerprintsTask = std::async(std::launch::async, distributeFingerprintTables, std::cref(args));
    size_t moleculesNumber = smilesTask.get();
    fingerprintsTask.get();
    timeTicker.tick("Files distribution");
    qtr::BallTreeBuilder ballTree(args.treeDepth, args.subtreeParallelDepth, args.dbDataDirsPaths,
                                  qtr::MaxDispersionBitSelector());
    timeTicker.tick("Ball tree building");
    std::ofstream ballTreeWriter(args.ballTreePath);
    ballTree.dumpNodes(ballTreeWriter);
    timeTicker.tick("Ball tree dumping");

    LOG(INFO) << "Molecules number: " << moleculesNumber;
    timeTicker.logResults();

    return 0;
}