#include "buildQtr.h"

#include "BallTreeBuilder.h"
#include "string_table_io/StringTableWriter.h"
#include "string_table_io/StringTableReader.h"
#include "fingerprint_table_io/FingerprintTableWriter.h"
#include "fingerprint_table_io/FingerprintTableReader.h"
#include "split_bit_selection/MaxDispersionBitSelector.h"
#include "properties_table_io/PropertiesTableReader.h"
#include "properties_table_io/PropertiesTableWriter.h"
#include "HuffmanCoder.h"
#include "PropertiesFilter.h"


#include <future>
#include <random>

using namespace std;

namespace qtr {

    namespace {
        void initFileSystem(const Args &args, TimeMeasurer &statisticCollector) {
            TimeMeasurer::FunctionExecutionTimer timer(statisticCollector, "filesystem initialization");

            vector<filesystem::path> alreadyExists;
            for (auto &dbDirPath: args.dbDataDirPaths()) {
                if (filesystem::exists(dbDirPath))
                    alreadyExists.emplace_back(dbDirPath);
            }
            if (filesystem::exists(args.dbOtherDataPath())) {
                alreadyExists.emplace_back(args.dbOtherDataPath());
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
            for (auto &dbDirPath: args.dbDataDirPaths()) {
                filesystem::create_directory(dbDirPath);
                filesystem::create_directory(dbDirPath / "0");
            }
            filesystem::create_directory(args.dbOtherDataPath());
            filesystem::create_directory(args.idToStringDestinationDirPath());
        }

        void distributeFingerprintTables(const Args &args, TimeMeasurer &statisticCollector) {
            TimeMeasurer::FunctionExecutionTimer timer(statisticCollector, "fingerprint tables distribution");

            vector<filesystem::path> ftFilePaths = findFiles(args.fingerprintTablesSourceDirPath(),
                                                             fingerprintTableExtension);
            shuffle(ftFilePaths.begin(), ftFilePaths.end(), mt19937(0));
            size_t drivesCount = args.dbDataDirPaths().size();
            for (size_t i = 0; i < ftFilePaths.size(); i++) {
                auto sourcePath = ftFilePaths[i];
                auto destinationPath = args.dbDataDirPaths()[i % drivesCount] / "0" / sourcePath.filename();
                filesystem::create_directory(destinationPath.parent_path());
                LOG(INFO) << "Copy " << sourcePath << " to " << destinationPath;
                filesystem::copy_file(sourcePath, destinationPath);
            }
        }

        map <uint64_t, filesystem::path> getLeafDirLocations(const Args &args) {
            map<uint64_t, filesystem::path> leafLocations;
            for (const auto &dirPath: args.dbDataDirPaths()) {
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

            auto dataDirPaths = args.dbDataDirPaths();
            for (size_t i = 0; i < leafIds.size(); i++) {
                const filesystem::path &currentLocation = leafLocations[leafIds[i]];
                const filesystem::path &newLocation = dataDirPaths[i % dataDirPaths.size()];
                filesystem::rename(currentLocation, newLocation);
            }
        }

        void buildBallTree(const Args &args, TimeMeasurer &statisticCollector) {
            TimeMeasurer::FunctionExecutionTimer timer(statisticCollector, "ball tree building");

            distributeFingerprintTables(args, statisticCollector);

            BallTreeBuilder ballTree(args.treeDepth(), args.parallelizeDepth(), args.dbDataDirPaths(),
                                     MaxDispersionBitSelector());
            ofstream ballTreeWriter(args.ballTreePath());
            ballTree.dumpNodes(ballTreeWriter);

            if (args.dbType() == Args::DataBaseType::QtrDrive)
                shuffleBallTreeLeaves(args);
        }

        void copyIdToStringTables(TimeMeasurer &statisticCollector, const filesystem::path &source,
                                  const filesystem::path &destination) {
            TimeMeasurer::FunctionExecutionTimer timer(statisticCollector,
                                                       "copying from " + source.string() + " to " +
                                                       destination.string());

            for (const filesystem::path &filePath: findFiles(source, ".csv")) {
                LOG(INFO) << "copy file" << filePath << " to " << destination / filePath.filename();
                filesystem::copy_file(filePath, destination / filePath.filename());
            }
        }

        size_t mergeSmilesTablesAndBuildHuffman(const Args &args, TimeMeasurer &statisticCollector) {
            TimeMeasurer::FunctionExecutionTimer timer(statisticCollector,
                                                       "merging smiles tables and building huffman");

            vector<filesystem::path> smilesTablePaths = findFiles(args.smilesSourceDirPath(), stringTableExtension);

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
            StringTableWriter writer(args.smilesTablePath());
            copy(smilesTable.begin(), smilesTable.end(), writer.begin());

            auto huffmanCoder = huffmanBuilder.build();
            huffmanCoder.dump(args.huffmanCoderPath());
            LOG(INFO) << "Finish merging smiles tables and building huffman";

            return smilesTable.size();
        }

        size_t mergePropertyTables(const Args &args, TimeMeasurer &statisticCollector) {
            TimeMeasurer::FunctionExecutionTimer timer(statisticCollector, "merging property tables");

            vector<filesystem::path> propertyTablePaths = findFiles(args.propertyTablesSourceDirPath(), "");

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
            PropertiesTableWriter writer(args.propertyTableDestinationPath());
            copy(propertyTable.begin(), propertyTable.end(), writer.begin());

            return propertyTable.size();
        }

        size_t mergeTables(const Args &args, TimeMeasurer &statisticCollector) {
            auto mergeSmilesTablesTask = async(launch::async, mergeSmilesTablesAndBuildHuffman, cref(args),
                                               ref(statisticCollector));
            auto mergePropertyTablesTask = async(launch::async, mergePropertyTables, cref(args),
                                                 ref(statisticCollector));

            size_t moleculesNumber = mergeSmilesTablesTask.get();
            size_t moleculesNumber2 = mergePropertyTablesTask.get();

            assert(moleculesNumber == moleculesNumber2);

            return moleculesNumber;
        }

        size_t distributeSmilesTables(const Args &args, const map <uint64_t, filesystem::path> &molLocations) {
            unordered_map<string, vector<pair<uint64_t, string>>> smilesTables;
            vector<filesystem::path> smilesTablePaths = findFiles(args.smilesSourceDirPath(), stringTableExtension);
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


        size_t distributePropertyTables(const Args &args, const map <uint64_t, filesystem::path> &molLocations) {
            unordered_map<string, vector<pair<uint64_t, PropertiesFilter::Properties>>> propertyTables;
            vector<filesystem::path> propertyTablePaths = findFiles(args.propertyTablesSourceDirPath(), "");
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

        size_t distributeTablesToLeafDirectories(const Args &args, TimeMeasurer &statisticCollector) {
            TimeMeasurer::FunctionExecutionTimer timer(statisticCollector, "smiles+properties tables distribution");

            auto leafLocations = getLeafDirLocations(args);

            map<uint64_t, filesystem::path> molLocations;
            for (auto &[leafId, leafDir]: leafLocations) {
                auto ftPath = leafDir / "data.ft";
                assert(filesystem::is_regular_file(ftPath));
                for (const auto &[molId, fp]: FingerprintTableReader(ftPath)) {
                    molLocations[molId] = leafDir;
                }
            }
            auto distributeSmilesTablesTask = async(launch::async, distributeSmilesTables, cref(args),
                                                    cref(molLocations));
            auto distributePropertyTablesTask = async(launch::async, distributePropertyTables, cref(args),
                                                      cref(molLocations));

            size_t moleculesNumber = distributeSmilesTablesTask.get();
            size_t moleculesNumber2 = distributePropertyTablesTask.get();
            assert(moleculesNumber == moleculesNumber2);
            return moleculesNumber;
        }
    }


    void buildQtrDb(const Args &args, TimeMeasurer &statisticCollector) {
        TimeMeasurer::FunctionExecutionTimer timer(statisticCollector, "db building");

        initFileSystem(args, statisticCollector);

        auto buildBallTreeTask = async(launch::async, buildBallTree, cref(args), ref(statisticCollector));

        auto idToStrSourceDirPath = args.idToStringSourceDirPath();
        auto idToStrDestinationDirPath = args.idToStringDestinationDirPath();
        auto copyIdToStrTablesTask = async(launch::async, copyIdToStringTables,
                                           ref(statisticCollector),
                                           cref(idToStrSourceDirPath),
                                           cref(idToStrDestinationDirPath));
        future<size_t> processTablesTask;
        if (args.dbType() == Args::DataBaseType::QtrRam)
            processTablesTask = async(launch::async, mergeTables, cref(args), ref(statisticCollector));
        else {
            buildBallTreeTask.wait(); // should distribute only after ball tree is built
            processTablesTask = async(launch::async, distributeTablesToLeafDirectories, cref(args),
                                      ref(statisticCollector));
        }

        buildBallTreeTask.wait();
        copyIdToStrTablesTask.wait();
        size_t moleculesNumber = processTablesTask.get();

        LOG(INFO) << "Molecules number: " << moleculesNumber;
    }
}