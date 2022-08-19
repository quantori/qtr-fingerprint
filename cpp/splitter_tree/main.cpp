#include "indigo.h"

#include "IndigoWriteBuffer.h"

#include <glog/logging.h>

#include <chrono>
#include <iostream>
#include <filesystem>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

#include "RawBucketsIO.h"
#include "CSVRawBucketIO.h"
#include "SplitterTree.h"
#include "ColumnsSelection.h"

ABSL_FLAG(std::vector<std::string>, rb_dir_paths, {},
          "Path to directory where raw bucket files to build structure are stored");

ABSL_FLAG(std::vector<std::string>, store_dir_paths, {},
          "Path to directory where data should be stored");

ABSL_FLAG(std::string, other_data_path, "",
          "Path to file where splitter tree should be stored");

ABSL_FLAG(std::string, cols_subset_path, "",
          "Path to file with columns for dimension reduction");

ABSL_FLAG(uint64_t, max_tree_depth, 0,
          "Max possible splitter tree depth");

ABSL_FLAG(uint64_t, stop_bucket_size, 0,
          "Node doesn't split if it has less than stop_bucket_size molecules inside");

ABSL_FLAG(uint64_t, subtree_parallel_depth, 0,
          "Depth on which subtree parallelization starts");

ABSL_FLAG(std::string, db_name, "",
          "Name of folders with data base's files");


struct Args {
    std::vector<std::filesystem::path> rbDirPaths;
    std::vector<std::filesystem::path> storeDirPaths;
    std::filesystem::path otherDataPath;
    std::string dbName;
    std::vector<std::filesystem::path> dbDataDirsPaths;
    std::filesystem::path columnsSubsetPath;
    uint64_t maxTreeDepth;
    uint64_t stopBucketSize;
    uint64_t subtreeParallelDepth;
    std::filesystem::path splitterTreePath;

    Args(int argc, char *argv[]) {
        absl::ParseCommandLine(argc, argv);

        std::vector<std::string> rbDirPathsStrings = absl::GetFlag(FLAGS_rb_dir_paths);
        std::copy(rbDirPathsStrings.begin(), rbDirPathsStrings.end(), std::back_inserter(rbDirPaths));
        qtr::emptyArgument(rbDirPaths, "Please specify rb_dir_paths option");
        for (size_t i = 0; i < rbDirPaths.size(); i++) {
            LOG(INFO) << "rbDirPaths[" << i << "]: " << rbDirPaths[i];
        }

        std::vector<std::string> storeDirPathsStrings = absl::GetFlag(FLAGS_store_dir_paths);
        std::copy(storeDirPathsStrings.begin(), storeDirPathsStrings.end(), std::back_inserter(storeDirPaths));
        qtr::emptyArgument(storeDirPaths, "Please specify store_dir_paths option");
        for (size_t i = 0; i < storeDirPaths.size(); i++) {
            LOG(INFO) << "storeDirPaths[" << i << "]: " << storeDirPaths[i];
        }

        otherDataPath = absl::GetFlag(FLAGS_other_data_path);
        qtr::emptyArgument(otherDataPath, "Please specify other_data_path option");
        LOG(INFO) << "otherDataPath: " << otherDataPath;

        dbName = absl::GetFlag(FLAGS_db_name);
        if (dbName.empty()) {
            LOG(INFO) << "db_name option is not specified";
            for (size_t i = 0; i < 10000; i++) {
                bool ok = true;
                std::string newDbName = "db_" + std::to_string(i);
                for (auto &dir: storeDirPaths) {
                    ok &= !std::filesystem::exists(dir / newDbName);
                }
                ok &= !std::filesystem::exists(otherDataPath / newDbName);
                if (ok) {
                    dbName = newDbName;
                    break;
                }
            }
            qtr::emptyArgument(dbName, "Something wrong with data base name generation");
            LOG(INFO) << "Generated name for data base: " << dbName;
        }
        LOG(INFO) << "dbName: " << dbName;

        for (auto &dir: storeDirPaths) {
            dbDataDirsPaths.emplace_back(dir / dbName);
        }
        for (size_t i = 0; i < dbDataDirsPaths.size(); i++) {
            LOG(INFO) << "dbDataDirPaths[" << i << "]: " << dbDataDirsPaths[i];
        }


        if (storeDirPaths.size() != rbDirPaths.size()) {
            LOG(ERROR) << "count of store dir paths have to be equal to count of rb dir paths";
            exit(-1);
        }

        splitterTreePath = otherDataPath / dbName / "tree";
        LOG(INFO) << "splitterTreePath: " << splitterTreePath;


        columnsSubsetPath = absl::GetFlag(FLAGS_cols_subset_path);
        qtr::emptyArgument(columnsSubsetPath, "Please specify cols_subset_path option");
        LOG(INFO) << "ColumnsSubsetPath: " << columnsSubsetPath;

        maxTreeDepth = absl::GetFlag(FLAGS_max_tree_depth);
        qtr::emptyArgument(maxTreeDepth, "Please specify max_tree_depth option");
        LOG(INFO) << "maxTreeDepth: " << maxTreeDepth;

        stopBucketSize = absl::GetFlag(FLAGS_stop_bucket_size);
        qtr::emptyArgument(stopBucketSize, "Please specify stop_bucket_size option");
        LOG(INFO) << "stopBucketSize: " << stopBucketSize;

        subtreeParallelDepth = absl::GetFlag(FLAGS_subtree_parallel_depth);
        qtr::emptyArgument(subtreeParallelDepth, "Please specify subtree_parallel_depth option");
        LOG(INFO) << "subtreeParallelDepth: " << subtreeParallelDepth;
    }
};

void askAboutContinue(const std::string &question) {
    std::cout << question << ". Continue? (Y/n): ";
    std::string userAnswer;
    std::cin >> userAnswer;
    if (userAnswer != "Y" && userAnswer != "y") {
        std::cout << "abort\n";
        exit(-1);
    }
}

void initFileSystem(const Args &args) {
    std::vector<std::filesystem::path> alreadyExists;
    for (auto &dbDirPath: args.dbDataDirsPaths) {
        if (!std::filesystem::create_directory(dbDirPath))
            alreadyExists.emplace_back(dbDirPath);
        std::filesystem::create_directory(dbDirPath / "0");
    }
    if (!alreadyExists.empty()) {
        std::cout << "Some store directories already exist: \n";
        for (auto &dir: alreadyExists) {
            std::cout << dir << '\n';
        }
        askAboutContinue("Data could be overridden");
    }
    if (!std::filesystem::create_directory(args.otherDataPath / args.dbName)) {
        std::cout << "Other file path " << args.otherDataPath / args.dbName << " already exists. ";
        askAboutContinue("Data could be overridden.");
    }


    assert(args.dbDataDirsPaths.size() == args.rbDirPaths.size());
    for (size_t i = 0; i < args.dbDataDirsPaths.size(); i++) {
        std::filesystem::create_directory(args.dbDataDirsPaths[i] / "0");
        for (auto &path: qtr::findFiles(args.rbDirPaths[i], ".rb")) {
            std::filesystem::copy(path, args.dbDataDirsPaths[i] / "0");
        }
    }
}

void initLogging(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    google::SetLogDestination(google::INFO, "info.log");
    FLAGS_alsologtostderr = true;
}

int main(int argc, char *argv[]) {

    initLogging(argc, argv);
    Args args(argc, argv);

    auto timePoints = std::vector{std::chrono::high_resolution_clock::now()};

    auto tickTimePoint = [&timePoints](const std::string &message) {
        auto t1 = timePoints.back();
        timePoints.emplace_back(std::chrono::high_resolution_clock::now());
        auto t2 = timePoints.back();
        std::chrono::duration<double> t = t2 - t1;
        LOG(INFO) << message << " : " << t.count() << " sec";
        return t.count();
    };

    initFileSystem(args);
    double initFileSystemTime = tickTimePoint("Init files");

    std::cerr << "db data dir paths: \n";
    for (auto &dir: args.dbDataDirsPaths) {
        std::cerr << dir << '\n';
    }

    qtr::SplitterTree tree(args.dbDataDirsPaths);
    tree.build(args.maxTreeDepth, args.stopBucketSize, args.subtreeParallelDepth);
    std::ofstream treeFileOut(args.otherDataPath);
    tree.dump(treeFileOut);
    double splitterTreeTime = tickTimePoint("Splitter tree is built");

    auto columnsSubset = qtr::ColumnsReader(args.columnsSubsetPath).readAll();
    auto selectFunction = qtr::PearsonCorrelationSelectionFunction(columnsSubset);
    auto columnsChooser = qtr::ColumnsSelector(args.dbDataDirsPaths, selectFunction);
    columnsChooser.handleRawBuckets();
    double columnsSelectingTime = tickTimePoint("Columns are selected");

    std::chrono::duration<double> elapsed_seconds = timePoints.back() - timePoints.front();
    LOG(INFO) << "Elapsed time: " << elapsed_seconds.count() << "s";
    LOG(INFO) << "Init files: " << initFileSystemTime << "s";
    LOG(INFO) << "Splitter tree building: " << splitterTreeTime << "s";
    LOG(INFO) << "Columns selecting: " << columnsSelectingTime << "s";
    return 0;
}
