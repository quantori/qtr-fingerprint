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
#include "HuffmanCoder.h"

using namespace std;
using namespace qtr;

ABSL_FLAG(string, smiles_dir_path, "",
          "Path to directory where smiles tables are stored");

ABSL_FLAG(string, fingerprints_dir_path, "",
          "Path to directory where fingerprint tables are stored");

ABSL_FLAG(vector<string>, data_dir_paths, {},
          "Path to directories where data should be stored");

ABSL_FLAG(string, other_data_path, "",
          "Path to file where splitter tree should be stored");

ABSL_FLAG(uint64_t, subtree_parallel_depth, 0,
          "Depth on which subtree parallelization starts");

ABSL_FLAG(uint64_t, tree_depth, 0,
          "ball tree depth");

ABSL_FLAG(string, db_name, "",
          "Name of folders with data base's files");

namespace {
    mt19937 randomGenerator(0);
}

struct Args {
    filesystem::path smilesDirPath;
    filesystem::path fingerprintsDirPath;
    vector<filesystem::path> dataDirPaths;
    filesystem::path otherDataPath;
    uint64_t subtreeParallelDepth;
    uint64_t treeDepth;
    string dbName;

    vector<filesystem::path> dbDataDirsPaths;
    filesystem::path dbOtherDataPath;

    filesystem::path ballTreePath;
    filesystem::path smilesTablePath;
    filesystem::path huffmanCoderPath;

    Args(int argc, char *argv[]) {
        absl::ParseCommandLine(argc, argv);

        smilesDirPath = absl::GetFlag(FLAGS_smiles_dir_path);
        emptyArgument(smilesDirPath, "Please specify smiles_dir_path option");
        LOG(INFO) << "smilesDirPath: " << smilesDirPath;

        fingerprintsDirPath = absl::GetFlag(FLAGS_fingerprints_dir_path);
        emptyArgument(fingerprintsDirPath, "Please specify fingerprints_dir_path option");
        LOG(INFO) << "fingerprintsDirPath" << fingerprintsDirPath;

        vector<string> dataDirPathsStrings = absl::GetFlag(FLAGS_data_dir_paths);
        copy(dataDirPathsStrings.begin(), dataDirPathsStrings.end(), back_inserter(dataDirPaths));
        emptyArgument(dataDirPaths, "Please specify data_dir_paths option");
        for (size_t i = 0; i < dataDirPaths.size(); i++) {
            LOG(INFO) << "dataDirPaths[" << i << "]: " << dataDirPaths[i];
        }

        otherDataPath = absl::GetFlag(FLAGS_other_data_path);
        emptyArgument(otherDataPath, "Please specify other_data_path option");
        LOG(INFO) << "otherDataPath: " << otherDataPath;

        subtreeParallelDepth = absl::GetFlag(FLAGS_subtree_parallel_depth);
        emptyArgument(subtreeParallelDepth, "Please specify subtree_parallel_depth option");
        LOG(INFO) << "subtreeParallelDepth: " << subtreeParallelDepth;

        treeDepth = absl::GetFlag(FLAGS_tree_depth);
        emptyArgument(treeDepth, "Please specify tree_depth option");
        LOG(INFO) << "treeDepth: " << treeDepth;

        dbName = absl::GetFlag(FLAGS_db_name);
        if (dbName.empty()) {
            LOG(INFO) << "db_name option is not specified";
            dbName = generateDbName(dataDirPaths, otherDataPath);
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

        huffmanCoderPath = dbOtherDataPath / "huffman";
        LOG(INFO) << "huffmanCoderPath: " << huffmanCoderPath;
    }
};

void initFileSystem(const Args &args) {
    vector<filesystem::path> alreadyExists;
    for (auto &dbDirPath: args.dbDataDirsPaths) {
        if (filesystem::exists(dbDirPath))
            alreadyExists.emplace_back(dbDirPath);
    }
    if (filesystem::exists(args.dbOtherDataPath)) {
        alreadyExists.emplace_back(args.dbOtherDataPath);
    }
    if (!alreadyExists.empty()) {
        cout << "Some data directories already exist: \n";
        for (auto &dir: alreadyExists) {
            cout << dir << '\n';
        }
        askAboutContinue("Data will be overridden");
        for (auto &dir: alreadyExists) {
            filesystem::remove_all(dir);
        }
    }
    for (auto &dbDirPath: args.dbDataDirsPaths) {
        filesystem::create_directory(dbDirPath);
        filesystem::create_directory(dbDirPath / "0");
    }
    filesystem::create_directory(args.dbOtherDataPath);
}

void distributeFingerprintTables(const Args &args) {
    vector<filesystem::path> ftFilePaths = findFiles(args.fingerprintsDirPath,
                                                     fingerprintTableExtension);
    shuffle(ftFilePaths.begin(), ftFilePaths.end(), randomGenerator);
    size_t drivesCount = args.dbDataDirsPaths.size();
    for (size_t i = 0; i < ftFilePaths.size(); i++) {
        auto sourcePath = ftFilePaths[i];
        auto destinationPath = args.dbDataDirsPaths[i % drivesCount] / "0" / sourcePath.filename();
        filesystem::create_directory(destinationPath.parent_path());
        LOG(INFO) << "Copy " << sourcePath << " to " << destinationPath;
        filesystem::copy_file(sourcePath, destinationPath);
    }
}

HuffmanCoder buildHuffman(const Args &args, const vector<filesystem::path> &smilesTablePaths) {
    HuffmanCoder::Builder huffmanBuilder;
    for (auto &stFile: smilesTablePaths) {
        for (const auto &[_, smiles]: SmilesTableReader(stFile)) {
            huffmanBuilder += smiles;
        }
    }
    HuffmanCoder huffmanCoder = huffmanBuilder.build();
    huffmanCoder.dump(args.huffmanCoderPath);
    return huffmanBuilder.build();
}

size_t mergeSmilesTables(const Args &args, const vector<filesystem::path> &smilesTablePaths) {
    vector<smiles_table_value_t> smilesTable;

    for (auto &stFile: smilesTablePaths) {
        SmilesTableReader reader(stFile);
        copy(reader.begin(), reader.end(), back_inserter(smilesTable));
    }

    sort(smilesTable.begin(), smilesTable.end(),
         [](const smiles_table_value_t &a, const smiles_table_value_t &b) {
             return a.first < b.first;
         });
    SmilesTableWriter writer(args.smilesTablePath);
    copy(smilesTable.begin(), smilesTable.end(), writer.begin());
    return smilesTable.size();
}

void buildDB(const Args &args) {
    TimeTicker timeTicker;

    initFileSystem(args);
    timeTicker.tick("Filesystem initialization");

    std::vector<std::filesystem::path> smilesTablePaths = qtr::findFiles(args.smilesDirPath, qtr::smilesTableExtension);

    auto mergeSmilesTask = async(launch::async, mergeSmilesTables, cref(args), cref(smilesTablePaths));
    auto fingerprintsTask = async(launch::async, distributeFingerprintTables, cref(args));
    fingerprintsTask.get();
    size_t moleculesNumber = mergeSmilesTask.get();

    auto huffmanTask = async(launch::async, buildHuffman, cref(args), cref(smilesTablePaths));

    timeTicker.tick("Fingerprint table files distribution + smiles merge + huffman build");
    BallTreeBuilder ballTree(args.treeDepth, args.subtreeParallelDepth, args.dbDataDirsPaths,
                             MaxDispersionBitSelector());

    timeTicker.tick("Ball tree building");
    ofstream ballTreeWriter(args.ballTreePath);

    ballTree.dumpNodes(ballTreeWriter);
    timeTicker.tick("Ball tree dumping");

    huffmanTask.get();

    LOG(INFO) << "Molecules number: " << moleculesNumber;

    timeTicker.logResults();
}

int main(int argc, char *argv[]) {
    initLogging(argv, google::INFO, "build_db.log", true);
    Args args(argc, argv);

    buildDB(args);

    return 0;
}