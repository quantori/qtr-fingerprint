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
#include "fingerprint_table_io/FingerprintTableReader.h"
#include "BallTreeBuilder.h"
#include "ball_tree/split_bit_selection/MaxDispersionBitSelector.h"
#include "HuffmanCoder.h"
#include "properties_table_io/PropertiesTableReader.h"
#include "properties_table_io/PropertiesTableWriter.h"
#include "PropertiesFilter.h"
#include "DbConfig.h"

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

ABSL_FLAG(string, db_type, "",
          "Possible types: on_drive, in_ram_smiles, in_ram_molecules");

TimeMeasurer statisticCollector;

struct Args {
    filesystem::path smilesSourceDirPath;
    filesystem::path fingerprintTablesSourceDirPath;
    filesystem::path idToStringSourceDirPath;
    filesystem::path propertyTablesSourceDirPath;

    vector<filesystem::path> destDirPaths;
    filesystem::path otherDataPath;
    uint64_t subtreeParallelDepth;
    uint64_t treeDepth;
    string dbName;

    DbType dbType;

    vector<filesystem::path> dbDataDirsPaths;
    filesystem::path dbOtherDataPath;

    filesystem::path ballTreePath;
    filesystem::path smilesTablePath;
    filesystem::path idToStringDestinationDirPath;
    filesystem::path huffmanCoderPath;
    filesystem::path propertyTableDestinationPath;


    Args(int argc, char *argv[]) {
        absl::ParseCommandLine(argc, argv);

        // get flags
        filesystem::path sourceDirPath = absl::GetFlag(FLAGS_source_dir_path);
        vector<string> tmpDestDirPaths = absl::GetFlag(FLAGS_dest_dir_paths);
        otherDataPath = absl::GetFlag(FLAGS_other_dest_dir_path);
        subtreeParallelDepth = absl::GetFlag(FLAGS_subtree_parallel_depth);
        treeDepth = absl::GetFlag(FLAGS_tree_depth);
        dbName = absl::GetFlag(FLAGS_db_name);
        string dbTypeStr = absl::GetFlag(FLAGS_db_type);

        // check empty flags
        checkEmptyArgument(sourceDirPath, "Please specify source_dir_path option");
        checkEmptyArgument(tmpDestDirPaths, "Please specify dest_dir_paths option");
        checkEmptyArgument(otherDataPath, "Please specify other_dest_dir_path option");
        checkEmptyArgument(subtreeParallelDepth, "Please specify subtree_parallel_depth option");
        checkEmptyArgument(treeDepth, "Please specify tree_depth option");
        checkEmptyArgument(dbTypeStr, "Please specify db_type option");

        // init values
        smilesSourceDirPath = sourceDirPath / "smilesTables";
        fingerprintTablesSourceDirPath = sourceDirPath / "fingerprintTables";
        idToStringSourceDirPath = sourceDirPath / "idToStringTables";
        propertyTablesSourceDirPath = sourceDirPath / "propertyTables";
        copy(tmpDestDirPaths.begin(), tmpDestDirPaths.end(), back_inserter(destDirPaths));
        if (dbName.empty()) {
            dbName = generateDbName(destDirPaths, otherDataPath);
        }
        for (auto &dir: destDirPaths) {
            dbDataDirsPaths.emplace_back(dir / dbName);
        }
        dbOtherDataPath = otherDataPath / dbName;
        ballTreePath = dbOtherDataPath / "tree";
        smilesTablePath = dbOtherDataPath / "moleculesTable";
        huffmanCoderPath = dbOtherDataPath / "huffman";
        idToStringDestinationDirPath = dbOtherDataPath / "idToString";
        propertyTableDestinationPath = dbOtherDataPath / "propertyTable";

        dbType = parseDbType(dbTypeStr);

        // log
        LOG(INFO) << "smilesSourceDirPath: " << smilesSourceDirPath;
        LOG(INFO) << "fingerprintTablesSourceDirPath" << fingerprintTablesSourceDirPath;
        LOG(INFO) << "idToStringSourceDirPath" << idToStringSourceDirPath;
        for (size_t i = 0; i < destDirPaths.size(); i++) {
            LOG(INFO) << "destDirPaths[" << i << "]: " << destDirPaths[i];
        }
        LOG(INFO) << "otherDataPath: " << otherDataPath;
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
    TimeMeasurer::FunctionTimeMeasurer timer(statisticCollector, "filesystem initialization");

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
    TimeMeasurer::FunctionTimeMeasurer timer(statisticCollector, "fingerprint tables distribution");

    vector<filesystem::path> ftFilePaths = findFiles(args.fingerprintTablesSourceDirPath,
                                                     fingerprintTableExtension);
    shuffle(ftFilePaths.begin(), ftFilePaths.end(), mt19937(0));
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
    TimeMeasurer::FunctionTimeMeasurer timer(statisticCollector, "merging smiles tables and building huffman");

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
    TimeMeasurer::FunctionTimeMeasurer timer(statisticCollector, "merging property tables");

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

    return propertyTable.size();
}

void copyIdToStringTables(const filesystem::path &source, const filesystem::path &destination) {
    TimeMeasurer::FunctionTimeMeasurer timer(statisticCollector,
                                             "copying from " + source.string() + " to " + destination.string());

    for (const filesystem::path &filePath: findFiles(source, ".csv")) {
        LOG(INFO) << "copy file" << filePath << " to " << destination / filePath.filename();
        filesystem::copy_file(filePath, destination / filePath.filename());
    }
}

size_t mergeTables(const Args &args) {
    auto mergeSmilesTablesTask = async(launch::async, mergeSmilesTablesAndBuildHuffman, cref(args));
    auto mergePropertyTablesTask = async(launch::async, mergePropertyTables, cref(args));

    size_t moleculesNumber = mergeSmilesTablesTask.get();
    size_t moleculesNumber2 = mergePropertyTablesTask.get();

    assert(moleculesNumber == moleculesNumber2);

    return moleculesNumber;
}

map<uint64_t, filesystem::path> getLeafDirLocations(const Args &args) {
    map<uint64_t, filesystem::path> leafLocations;
    for (const auto &dirPath: args.dbDataDirsPaths) {
        for (auto &filePath: findFiles(dirPath)) {
            uint64_t leafId = stoi(filePath.stem());
            leafLocations[leafId] = filePath;
        }
    }
    return leafLocations;
}

void shuffleBallTreeLeaves(const Args &args) {
    auto leafLocations = getLeafDirLocations(args);

    vector<uint64_t> leafIds;
    for (auto &[id, _]: leafLocations)
        leafIds.push_back(id);
    shuffle(leafIds.begin(), leafIds.end(), mt19937(0));

    for (size_t i = 0; i < leafIds.size(); i++) {
        const filesystem::path &currentLocation = leafLocations[leafIds[i]];
        const filesystem::path &newLocation = args.dbDataDirsPaths[i % args.dbDataDirsPaths.size()];
        filesystem::rename(currentLocation, newLocation);
    }
}

size_t distributeSmilesTables(const Args &args, const map<uint64_t, filesystem::path> &molLocations) {
    unordered_map<string, vector<pair<uint64_t, string>>> smilesTables;
    vector<filesystem::path> smilesTablePaths = findFiles(args.smilesSourceDirPath, stringTableExtension);
    size_t molNumber = 0;
    for (auto &tablePath: smilesTablePaths) {
        for (const auto &[id, smiles]: StringTableReader(tablePath)) {
            const auto &location = molLocations.find(id)->second;
            smilesTables[location].emplace_back(id, smiles);
            molNumber++;
        }
    }

    for (auto &[leafDirPath, leafTable]: smilesTables) {
        auto leafTablePath = filesystem::path(leafDirPath) / ("smiles" + stringTableExtension);
        StringTableWriter writer(leafTablePath);
        writer << leafTable;
    }

    return molNumber;
}


size_t distributePropertyTables(const Args &args, const map<uint64_t, filesystem::path> &molLocations) {
    unordered_map<string, vector<pair<uint64_t, PropertiesFilter::Properties>>> propertyTables;
    vector<filesystem::path> propertyTablePaths = findFiles(args.propertyTablesSourceDirPath, "");
    size_t molNumber = 0;
    for (auto &tablePath: propertyTablePaths) {
        for (const auto &[id, properties]: PropertiesTableReader(tablePath)) {
            const auto &location = molLocations.find(id)->second;
            propertyTables[location].emplace_back(id, properties);
            molNumber++;
        }
    }

    for (auto &[leafDirStr, leafTable]: propertyTables) {
        auto leafTablePath = filesystem::path(leafDirStr) / "properties";
        PropertiesTableWriter writer(leafTablePath);
        writer << leafTable;
    }

    return molNumber;
}

size_t distributeTablesToLeafDirectories(const Args &args) {
    TimeMeasurer::FunctionTimeMeasurer timer(statisticCollector, "smiles+properties tables distribution");

    auto leafLocations = getLeafDirLocations(args);

    map<uint64_t, filesystem::path> molLocations;
    for (auto &[leafId, leafDir]: leafLocations) {
        auto ftPath = leafDir / "data.ft";
        assert(filesystem::is_regular_file(ftPath));
        for (const auto &[molId, fp]: FingerprintTableReader(ftPath)) {
            molLocations[molId] = leafDir;
        }
    }
    auto distributeSmilesTablesTask = async(launch::async, distributeSmilesTables, cref(args), cref(molLocations));
    auto distributePropertyTablesTask = async(launch::async, distributePropertyTables, cref(args), cref(molLocations));

    size_t moleculesNumber = distributeSmilesTablesTask.get();
    size_t moleculesNumber2 = distributePropertyTablesTask.get();
    assert(moleculesNumber == moleculesNumber2);
    return moleculesNumber;
}

void buildBallTree(const Args &args) {
    TimeMeasurer::FunctionTimeMeasurer timer(statisticCollector, "ball tree building");

    distributeFingerprintTables(args);

    BallTreeBuilder ballTree(args.treeDepth, args.subtreeParallelDepth, args.dbDataDirsPaths,
                             MaxDispersionBitSelector());
    ofstream ballTreeWriter(args.ballTreePath);
    ballTree.dumpNodes(ballTreeWriter);

    if (args.dbType == DbType::OnDrive)
        shuffleBallTreeLeaves(args);
}

void buildDb(const Args &args) {
    TimeMeasurer::FunctionTimeMeasurer timer(statisticCollector, "db building");

    initFileSystem(args);

    auto buildBallTreeTask = async(launch::async, buildBallTree, cref(args));
    auto copyIdToStrTablesTask = async(launch::async, copyIdToStringTables,
                                       cref(args.idToStringSourceDirPath),
                                       cref(args.idToStringDestinationDirPath));
    future<size_t> processTablesTask;
    if (args.dbType == DbType::InRamSmiles || args.dbType == DbType::InRamMolecules)
        processTablesTask = async(launch::async, mergeTables, cref(args));
    else {
        buildBallTreeTask.wait(); // should distribute only after ball tree is built
        processTablesTask = async(launch::async, distributeTablesToLeafDirectories, cref(args));
    }

    buildBallTreeTask.wait();
    copyIdToStrTablesTask.wait();
    size_t moleculesNumber = processTablesTask.get();

    LOG(INFO) << "Molecules number: " << moleculesNumber;
}

int main(int argc, char *argv[]) {
    initLogging(argv, google::INFO, "build_db.info", true);
    Args args(argc, argv);

    buildDb(args);

    for (auto &[label, time]: statisticCollector) {
        LOG(INFO) << label << ": " << time;
    }

    return 0;
}