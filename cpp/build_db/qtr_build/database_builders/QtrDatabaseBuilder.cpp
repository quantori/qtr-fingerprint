#include "QtrDatabaseBuilder.h"

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
#include "Profiling.h"


#include <future>
#include <random>

using namespace std;

namespace qtr {
    namespace {
        void initFileSystem(const BuildArgs &args) {
            ProfileScope("Filesystem initialization");

            vector<filesystem::path> alreadyExists;
            for (auto &dbDir: args.dbDataDirs()) {
                if (filesystem::exists(dbDir))
                    alreadyExists.emplace_back(dbDir);
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
            for (auto &dbDir: args.dbDataDirs()) {
                filesystem::create_directory(dbDir);
                filesystem::create_directory(dbDir / "0");
            }
            filesystem::create_directory(args.dbOtherDataPath());
            filesystem::create_directory(args.idToStringDestinationDir());
        }

        void distributeFingerprintTables(const BuildArgs &args) {
            ProfileScope("Fingerprint tables distribution");

            vector<filesystem::path> ftFilePaths = findFiles(args.fingerprintTablesSourceDir(),
                                                             fingerprintTableExtension);
            shuffle(ftFilePaths.begin(), ftFilePaths.end(), mt19937(0));
            auto dbDataDirs = args.dbDataDirs();
            size_t drivesCount = dbDataDirs.size();
            for (size_t i = 0; i < ftFilePaths.size(); i++) {
                auto sourcePath = ftFilePaths[i];
                auto destinationPath = dbDataDirs[i % drivesCount] / "0" / sourcePath.filename();
                filesystem::create_directory(destinationPath.parent_path());
                LOG(INFO) << "Copy " << sourcePath << " to " << destinationPath;
                copyFileAndCheck(sourcePath, destinationPath);
            }
        }

        map <uint64_t, filesystem::path> getLeafDirLocations(const BuildArgs &args) {
            map<uint64_t, filesystem::path> leafLocations;
            for (const auto &dir: args.dbDataDirs()) {
                for (auto &filePath: findFiles(dir, "")) {
                    uint64_t leafId = stoi(filePath.stem());
                    leafLocations[leafId] = filePath;
                }
            }
            return leafLocations;
        }

        void shuffleBallTreeLeaves(const BuildArgs &args) {
            auto leafLocations = getLeafDirLocations(args);

            vector<uint64_t> leafIds;
            for (auto &[id, _]: leafLocations)
                leafIds.push_back(id);
            shuffle(leafIds.begin(), leafIds.end(), mt19937(0));

            auto dataDirs = args.dbDataDirs();
            for (size_t i = 0; i < leafIds.size(); i++) {
                const filesystem::path &currentLocation = leafLocations[leafIds[i]];
                const filesystem::path &newLocation = dataDirs[i % dataDirs.size()];
                filesystem::rename(currentLocation, newLocation);
            }
        }

        void buildBallTree(const BuildArgs &args) {
            ProfileScope("Ball tree building");

            distributeFingerprintTables(args);

            BallTreeBuilder ballTree(args.treeDepth(), args.parallelizeDepth(), args.dbDataDirs(),
                                     MaxDispersionBitSelector());
            ofstream ballTreeWriter(args.ballTreePath());
            ballTree.dumpNodes(ballTreeWriter);

            if (args.dbType() == DatabaseType::QtrDrive)
                shuffleBallTreeLeaves(args);
        }

        void copyIdToStringTables(const filesystem::path &source, const filesystem::path &destination) {
            ProfileScope("copying from " + source.string() + " to " + destination.string());

            for (const filesystem::path &filePath: findFiles(source, ".csv")) {
                LOG(INFO) << "copy file" << filePath << " to " << destination / filePath.filename();
                copyFileAndCheck(filePath, destination / filePath.filename());
            }
        }

        size_t mergeSmilesTablesAndBuildHuffman(const BuildArgs &args) {
            ProfileScope("Merging smiles tables and building huffman");

            vector<filesystem::path> smilesTablePaths = findFiles(args.smilesSourceDir(), stringTableExtension);

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

        size_t mergePropertyTables(const BuildArgs &args) {
            ProfileScope("Merge property tables");

            vector<filesystem::path> propertyTablePaths = findFiles(args.propertyTablesSourceDir(), "");

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

        size_t mergeTables(const BuildArgs &args) {
            auto mergeSmilesTablesTask = async(launch::async, mergeSmilesTablesAndBuildHuffman, cref(args));

            if (args.properties()) {
                auto mergePropertyTablesTask = async(launch::async, mergePropertyTables, cref(args));
                mergePropertyTablesTask.wait();
            }
            size_t moleculesNumber = mergeSmilesTablesTask.get();
            return moleculesNumber;
        }

        size_t distributeSmilesTables(const BuildArgs &args, const map <uint64_t, filesystem::path> &molLocations) {
            unordered_map<string, vector<pair<uint64_t, string>>> smilesTables;
            vector<filesystem::path> smilesTablePaths = findFiles(args.smilesSourceDir(), stringTableExtension);
            size_t molNumber = 0;
            for (auto &tablePath: smilesTablePaths) {
                for (const auto &[id, smiles]: StringTableReader(tablePath)) {
                    const auto &location = molLocations.find(id)->second;
                    smilesTables[location].emplace_back(id, smiles);
                    molNumber++;
                }
            }

            for (auto &[leafDir, leafTable]: smilesTables) {
                auto leafTablePath = filesystem::path(leafDir) / ("smiles" + stringTableExtension);
                StringTableWriter writer(leafTablePath);
                writer << leafTable;
            }

            return molNumber;
        }


        size_t distributePropertyTables(const BuildArgs &args, const map <uint64_t, filesystem::path> &molLocations) {
            unordered_map<string, vector<pair<uint64_t, PropertiesFilter::Properties>>> propertyTables;
            vector<filesystem::path> propertyTablePaths = findFiles(args.propertyTablesSourceDir(), "");
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

        size_t distributeTablesToLeafDirectories(const BuildArgs &args) {
            ProfileScope("SMILES+Properties tables distribution");

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

            if (args.properties()) {
                auto distributePropertyTablesTask = async(launch::async, distributePropertyTables, cref(args),
                                                          cref(molLocations));
                distributePropertyTablesTask.wait();
            }
            size_t moleculesNumber = distributeSmilesTablesTask.get();
            return moleculesNumber;
        }
    }

    void QtrDatabaseBuilder::build(const BuildArgs &args) {
        ProfileScope("QtrDB building");

        initFileSystem(args);

        auto buildBallTreeTask = async(launch::async, buildBallTree, cref(args));

        auto idToStrSourceDir = args.idToStringSourceDir();
        auto idToStrDestinationDir = args.idToStringDestinationDir();
        auto copyIdToStrTablesTask = async(launch::async, copyIdToStringTables,
                                           cref(idToStrSourceDir),
                                           cref(idToStrDestinationDir));
        future<size_t> processTablesTask;
        if (args.dbType() == DatabaseType::QtrRam)
            processTablesTask = async(launch::async, mergeTables, cref(args));
        else if (args.dbType() == DatabaseType::QtrDrive) {
            buildBallTreeTask.wait(); // should distribute only after ball tree is built
            processTablesTask = async(launch::async, distributeTablesToLeafDirectories, cref(args));
        } else {
            LOG(ERROR) << "Wrong DB type in QtrDB building function";
            exit(-1);
        }

        buildBallTreeTask.wait();
        copyIdToStrTablesTask.wait();
        size_t moleculesNumber = processTablesTask.get();

        LOG(INFO) << "Molecules number: " << moleculesNumber;
    }
} // qtr