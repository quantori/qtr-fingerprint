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

using namespace indigo_cpp;
using namespace qtr;
using namespace std;

#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,fma,tune=native")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")

ABSL_FLAG(std::string, path_to_store_dir, "", "Path to directory where data should be stored");

ABSL_FLAG(std::string, path_to_rb_files_dir, "",
          "Path to directory where raw bucket files to build structure are stored");

ABSL_FLAG(std::string, path_to_cols_subset, "",
          "Path to file with columns for dimension reduction");

filesystem::path getStoreDirPath(const filesystem::path &dataDir) {
    for (size_t i = 0;; i++) {
        auto path = dataDir / ("search_data_" + to_string(i));
        if (filesystem::exists(path))
            continue;
        filesystem::create_directory(path);
        LOG(INFO) << "Store data to " << path << endl;
        return path;
    }
}

void initLogging(int argc, char *argv[]) {
    google::InitGoogleLogging(argv[0]);
    google::SetLogDestination(google::INFO, "info.log");
    FLAGS_alsologtostderr = true;
}

int main(int argc, char *argv[]) {

    initLogging(argc, argv);
    absl::ParseCommandLine(argc, argv);

    std::filesystem::path baseStoreDirPath = absl::GetFlag(FLAGS_path_to_store_dir);
    std::filesystem::path rbFilesPath = absl::GetFlag(FLAGS_path_to_rb_files_dir);
    std::filesystem::path columnsSubsetPath = absl::GetFlag(FLAGS_path_to_cols_subset);

    emptyArgument(baseStoreDirPath, "Please specify path_to_store_dir option");
    emptyArgument(rbFilesPath, "Please specify path_to_rb_files_dir option");
    if (columnsSubsetPath.empty()) {
        LOG(WARNING) << "Probably, you forgot to specify path_to_cols_subset option. Full list of columns will be used";
    }

    auto storeDirPath = getStoreDirPath(baseStoreDirPath);
    auto rawBucketsDirPath = storeDirPath / "raw_buckets";
    auto splitterTreeFilePath = storeDirPath / "tree";
    auto treeRootNodeDirPath = rawBucketsDirPath / "0";

    filesystem::create_directory(rawBucketsDirPath);
    filesystem::create_directory(treeRootNodeDirPath);
    filesystem::copy(rbFilesPath, treeRootNodeDirPath);

    auto timePoints = vector{chrono::high_resolution_clock::now()};

    auto tickTimePoint = [&timePoints](const string &message) {
        auto t1 = timePoints.back();
        timePoints.emplace_back(chrono::high_resolution_clock::now());
        auto t2 = timePoints.back();
        chrono::duration<double> t = t2 - t1;
        LOG(INFO) << message << " : " << t.count() << " sec" << endl;
    };

    // Build splitter tree
    SplitterTree tree(rawBucketsDirPath);
    auto bucketDirPaths = tree.build(13, 1000, 3); // 16, 600, 3
    ofstream treeFileOut(splitterTreeFilePath);
    tree.dump(treeFileOut);
    tickTimePoint("Splitter tree is built");

    // Choose minimum correlated columns
    auto columnsSubset = ColumnsReader(columnsSubsetPath).readAll();
    auto selectFunction = qtr::PearsonCorrelationSelectionFunction(columnsSubset);
    auto columnsChooser = ColumnsSelector(rawBucketsDirPath, selectFunction);
    columnsChooser.handleRawBuckets();
    tickTimePoint("Columns are chosen");
#pragma omp parallel for
    for (const auto &dirPath: bucketDirPaths) {
        assert(filesystem::is_directory(dirPath));
        LOG(INFO) << "Start creating CSV for " << dirPath;
        for (const auto &filePath: findFiles(dirPath, rawBucketExtension)) {
            {
                RawBucketReader reader(filePath);
                filesystem::path outFilePath = filePath;
                outFilePath.replace_extension(csvFileExtension);
                CSVRawBucketWriter writer(outFilePath);
                std::copy(reader.begin(), reader.end(), writer.begin());
            }
            filesystem::remove(filePath);
        }
        LOG(INFO) << "Finish creating CSV for " << dirPath;
    }

    chrono::duration<double> elapsed_seconds = timePoints.back() - timePoints.front();
    LOG(INFO) << "Elapsed time: " << elapsed_seconds.count() << "s\n";
    return 0;
}
