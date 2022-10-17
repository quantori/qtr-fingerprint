#include "indigo.h"

#include "IndigoWriteBuffer.h"

#include <glog/logging.h>

#include <chrono>
#include <iostream>
#include <filesystem>
#include <random>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

#include "raw_bucket_io/RawBucketReader.h"
#include "raw_bucket_io/RawBucketWriter.h"
#include "SplitterTree.h"
#include "ColumnsSelection.h"

ABSL_FLAG(std::string, rb_dir_path, {},
          "Path to directory where raw bucket files to build structure are stored");

ABSL_FLAG(std::vector<std::string>, data_dir_paths, {},
          "Path to directories where data should be stored");

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

ABSL_FLAG(std::string, raw_db_name, "",
          "Name of folders with data base's files");

namespace {
    std::mt19937 randomGenerator(0);
}


struct Args {
    std::filesystem::path rbDirPath;
    std::vector<std::filesystem::path> dataDirPaths;
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

        rbDirPath = absl::GetFlag(FLAGS_rb_dir_path);
        qtr::emptyArgument(rbDirPath, "Please specify rb_dir_path option");
        LOG(INFO) << "rbDirPath: " << rbDirPath;

        std::vector<std::string> dataDirPathsStrings = absl::GetFlag(FLAGS_data_dir_paths);
        std::copy(dataDirPathsStrings.begin(), dataDirPathsStrings.end(), std::back_inserter(dataDirPaths));
        qtr::emptyArgument(dataDirPaths, "Please specify data_dir_paths option");
        for (size_t i = 0; i < dataDirPaths.size(); i++) {
            LOG(INFO) << "dataDirPaths[" << i << "]: " << dataDirPaths[i];
        }

        otherDataPath = absl::GetFlag(FLAGS_other_data_path);
        qtr::emptyArgument(otherDataPath, "Please specify other_data_path option");
        LOG(INFO) << "otherDataPath: " << otherDataPath;

        dbName = absl::GetFlag(FLAGS_raw_db_name);
        if (dbName.empty()) {
            LOG(INFO) << "raw_db_name option is not specified";
            dbName = qtr::generateDbName(dataDirPaths, otherDataPath);
            LOG(INFO) << "Generated name for raw data base: " << dbName;
        }
        LOG(INFO) << "dbName: " << dbName;

        for (auto &dir: dataDirPaths) {
            dbDataDirsPaths.emplace_back(dir / dbName);
        }
        for (size_t i = 0; i < dbDataDirsPaths.size(); i++) {
            LOG(INFO) << "dbDataDirPaths[" << i << "]: " << dbDataDirsPaths[i];
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

void initFileSystem(const Args &args) {
    std::vector<std::filesystem::path> alreadyExists;
    for (auto &dbDirPath: args.dbDataDirsPaths) {
        if (std::filesystem::exists(dbDirPath))
            alreadyExists.emplace_back(dbDirPath);
    }
    if (std::filesystem::exists(args.otherDataPath / args.dbName)) {
        alreadyExists.emplace_back(args.otherDataPath / args.dbName);
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
    std::filesystem::create_directory(args.otherDataPath / args.dbName);

    std::vector<std::filesystem::path> rbFilePaths = qtr::findFiles(args.rbDirPath, ".rb");
    std::shuffle(rbFilePaths.begin(), rbFilePaths.end(), randomGenerator);
    size_t drivesCount = args.dbDataDirsPaths.size();
    for (size_t dataDirId = 0, rbFileId = 0; dataDirId < drivesCount; dataDirId++) {
        size_t currDriveFilesCount = rbFilePaths.size() / drivesCount
                                     + size_t(dataDirId < rbFilePaths.size() % drivesCount);
        for (size_t i = 0; i < currDriveFilesCount; i++, rbFileId++) {
            auto sourcePath = rbFilePaths[rbFileId];
            auto destinationPath = args.dbDataDirsPaths[dataDirId] / "0" / sourcePath.filename();
            std::filesystem::create_directory(destinationPath.parent_path());
            LOG(INFO) << "Copy " << sourcePath << " to " << destinationPath;
            std::filesystem::copy_file(sourcePath, destinationPath);
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

    qtr::TimeTicker timeTicker;

    initFileSystem(args);
    timeTicker.tick("Files initialization");

    qtr::SplitterTree tree(args.dbDataDirsPaths);
    tree.build(args.maxTreeDepth, args.stopBucketSize, args.subtreeParallelDepth);
    std::ofstream treeFileOut(args.splitterTreePath);
    tree.dump(treeFileOut);
    timeTicker.tick("Splitter tree building");

    std::vector<size_t> columnsSubset;
    qtr::IndexesReader(args.columnsSubsetPath) >> columnsSubset;
    auto selectFunction = qtr::PearsonCorrelationSelectionFunction(columnsSubset);
    auto columnsSelector = qtr::ColumnsSelector(args.dbDataDirsPaths, selectFunction);
    columnsSelector.handleRawBuckets();
    timeTicker.tick("Columns selection");

    timeTicker.logResults();
    return 0;
}
