#include <filesystem>
#include <string>
#include <random>
#include <future>

#include <glog/logging.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>

#include "Utils.h"
#include "string_table_io/StringTableWriter.h"
#include "string_table_io/StringTableReader.h"
#include "fingerprint_table_io/FingerprintTableWriter.h"
#include "BallTreeBuilder.h"
#include "ball_tree/split_bit_selection/MaxDispersionBitSelector.h"
#include "HuffmanCoder.h"
#include "properties_table_io/PropertiesTableReader.h"
#include "properties_table_io/PropertiesTableWriter.h"
#include "PropertiesFilter.h"

using namespace std;
using namespace qtr;

ABSL_FLAG(string, source_dir_path, "",
          "Path to directory where source data are stored");

ABSL_FLAG(vector<string>, dest_dir_paths, {},
          "Path to directories where data should be stored");

ABSL_FLAG(string, other_dest_dir_path, "",
          "Path to directory where other data should be stored");

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
    filesystem::path smilesSourceDirPath;
    filesystem::path fingerprintTablesSourceDirPath;
    filesystem::path idToStringSourceDirPath;
    filesystem::path propertyTablesSourceDirPath;

    vector<filesystem::path> destDirPaths;
    filesystem::path otherDestDirPath;
    uint64_t subtreeParallelDepth;
    uint64_t treeDepth;
    string dbName;

    vector<filesystem::path> dbDataDirsPaths;
    filesystem::path dbOtherDataPath;

    filesystem::path ballTreePath;
    filesystem::path smilesTablePath;
    filesystem::path idToStringDestinationDirPath;
    filesystem::path huffmanCoderPath;
    filesystem::path propertyTableDestinationPath;


    Args(int argc, char *argv[]) {
        absl::ParseCommandLine(argc, argv);

        filesystem::path sourceDirPath = absl::GetFlag(FLAGS_source_dir_path);
        vector<string> tmpDestDirPaths = absl::GetFlag(FLAGS_dest_dir_paths);
        otherDestDirPath = absl::GetFlag(FLAGS_other_dest_dir_path);
        subtreeParallelDepth = absl::GetFlag(FLAGS_subtree_parallel_depth);
        treeDepth = absl::GetFlag(FLAGS_tree_depth);
        dbName = absl::GetFlag(FLAGS_db_name);

        checkEmptyArgument(sourceDirPath, "Please specify source_dir_path option");
        checkEmptyArgument(tmpDestDirPaths, "Please specify dest_dir_paths option");
        checkEmptyArgument(otherDestDirPath, "Please specify other_dest_dir_path option");
        checkEmptyArgument(subtreeParallelDepth, "Please specify subtree_parallel_depth option");
        checkEmptyArgument(treeDepth, "Please specify tree_depth option");

        smilesSourceDirPath = sourceDirPath / "smilesTables";
        fingerprintTablesSourceDirPath = sourceDirPath / "fingerprintTables";
        idToStringSourceDirPath = sourceDirPath / "idToStringTables";
        propertyTablesSourceDirPath = sourceDirPath / "propertyTables";
        copy(tmpDestDirPaths.begin(), tmpDestDirPaths.end(), back_inserter(destDirPaths));
        if (dbName.empty()) {
            dbName = generateDbName(destDirPaths, otherDestDirPath);
        }
        for (auto &dir: destDirPaths) {
            dbDataDirsPaths.emplace_back(dir / dbName);
        }
        dbOtherDataPath = otherDestDirPath / dbName;
        ballTreePath = dbOtherDataPath / "tree";
        smilesTablePath = dbOtherDataPath / "smilesTable";
        huffmanCoderPath = dbOtherDataPath / "huffman";
        idToStringDestinationDirPath = dbOtherDataPath / "idToString";
        propertyTableDestinationPath = dbOtherDataPath / "propertyTable";

        LOG(INFO) << "smilesSourceDirPath: " << smilesSourceDirPath;
        LOG(INFO) << "fingerprintTablesSourceDirPath" << fingerprintTablesSourceDirPath;
        LOG(INFO) << "idToStringSourceDirPath" << idToStringSourceDirPath;
        for (size_t i = 0; i < destDirPaths.size(); i++) {
            LOG(INFO) << "destDirPaths[" << i << "]: " << destDirPaths[i];
        }
        LOG(INFO) << "otherDestDirPath: " << otherDestDirPath;
        LOG(INFO) << "subtreeParallelDepth: " << subtreeParallelDepth;
        LOG(INFO) << "treeDepth: " << treeDepth;
        LOG(INFO) << "dbName: " << dbName;
        for (size_t i = 0; i < dbDataDirsPaths.size(); i++) {
            LOG(INFO) << "dbDataDirPaths[" << i << "]: " << dbDataDirsPaths[i];
        }
        LOG(INFO) << "dbOtherDataPath: " << dbOtherDataPath;
        LOG(INFO) << "splitterTreePath: " << ballTreePath;
        LOG(INFO) << "smilesTablePath: " << smilesTablePath;
        LOG(INFO) << "huffmanCoderPath: " << huffmanCoderPath;
        LOG(INFO) << "idToStringDestinationDirPath: " << idToStringDestinationDirPath;
        LOG(INFO) << "propertyTableDestinationPath: " << propertyTableDestinationPath;

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
    filesystem::create_directory(args.idToStringDestinationDirPath);
}

void distributeFingerprintTables(const Args &args) {
    vector<filesystem::path> ftFilePaths = findFiles(args.fingerprintTablesSourceDirPath,
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

size_t mergeSmilesTablesAndBuildHuffman(const Args &args) {
    LOG(INFO) << "Start merging smiles tables and building huffman";
    vector<filesystem::path> smilesTablePaths = findFiles(args.smilesSourceDirPath, stringTableExtension);

    vector<string_table_value_t> smilesTable;
    HuffmanCoder::Builder huffmanBuilder;

    for (auto &stFile: smilesTablePaths) {
        StringTableReader reader(stFile);
        for (const auto &value: reader) {
            smilesTable.emplace_back(value);
            huffmanBuilder += value.second;
        }
    }

    sort(smilesTable.begin(), smilesTable.end(),
         [](const string_table_value_t &a, const string_table_value_t &b) {
             return a.first < b.first;
         });
    StringTableWriter writer(args.smilesTablePath);
    copy(smilesTable.begin(), smilesTable.end(), writer.begin());

    auto huffmanCoder = huffmanBuilder.build();
    huffmanCoder.dump(args.huffmanCoderPath);
    LOG(INFO) << "Finish merging smiles tables and building huffman";

    return smilesTable.size();
}

size_t mergePropertyTables(const Args &args) {
    LOG(INFO) << "Start merging property tables";
    vector<filesystem::path> propertyTablePaths = findFiles(args.propertyTablesSourceDirPath, "");

    vector<pair<uint64_t, PropertiesFilter::Properties>> propertyTable;

    for (auto &file: propertyTablePaths) {
        PropertiesTableReader reader(file);
        for (const auto &value: reader) {
            propertyTable.emplace_back(value);
        }
    }

    sort(propertyTable.begin(), propertyTable.end(),
         [](const pair<uint64_t, PropertiesFilter::Properties> &a,
            const pair<uint64_t, PropertiesFilter::Properties> &b) {
             return a.first < b.first;
         });
    for (size_t i = 0; i < propertyTable.size(); i++) {
        assert(propertyTable[i].first == i);
    }
    PropertiesTableWriter writer(args.propertyTableDestinationPath);
    copy(propertyTable.begin(), propertyTable.end(), writer.begin());
    LOG(INFO) << "Finish merging property tables";

    return propertyTable.size();
}

void copyIdToStringTables(const filesystem::path &source, const filesystem::path &destination) {
    LOG(INFO) << "Start copying from " << source << " to " << destination;
    for (const filesystem::path &filePath: findFiles(source, ".csv")) {
        LOG(INFO) << "copy file" << filePath << " to " << destination / filePath.filename();
        filesystem::copy_file(filePath, destination / filePath.filename());
    }
    LOG(INFO) << "Finish copying from " << source << " to " << destination;
}

size_t processTables(const Args &args) {
    auto copyIdToStrTablesTask = std::async(std::launch::async, copyIdToStringTables,
                                            cref(args.idToStringSourceDirPath),
                                            cref(args.idToStringDestinationDirPath));
    auto mergeSmilesTablesTask = async(launch::async, mergeSmilesTablesAndBuildHuffman, cref(args));
    auto fingerprintsTask = async(launch::async, distributeFingerprintTables, cref(args));
    auto mergePropertyTablesTask = async(launch::async, mergePropertyTables, cref(args));

    size_t moleculesNumber = mergeSmilesTablesTask.get();
    copyIdToStrTablesTask.get();
    fingerprintsTask.get();
    size_t moleculesNumber2 = mergePropertyTablesTask.get();

    assert(moleculesNumber == moleculesNumber2);

    return moleculesNumber;
}

void buildDB(const Args &args) {
    TimeTicker timeTicker;

    initFileSystem(args);
    timeTicker.tick("Filesystem initialization");
    size_t moleculesNumber = processTables(args);
    timeTicker.tick("Fingerprint table files distribution + smiles merge + huffman build");
    BallTreeBuilder ballTree(args.treeDepth, args.subtreeParallelDepth, args.dbDataDirsPaths,
                             MaxDispersionBitSelector());
    timeTicker.tick("Ball tree building");
    ofstream ballTreeWriter(args.ballTreePath);
    ballTree.dumpNodes(ballTreeWriter);
    timeTicker.tick("Ball tree dumping");

    LOG(INFO) << "Molecules number: " << moleculesNumber;

    timeTicker.logResults();
}

int main(int argc, char *argv[]) {
    initLogging(argv, google::INFO, "build_db.log", true);
    Args args(argc, argv);

    buildDB(args);

    return 0;
}