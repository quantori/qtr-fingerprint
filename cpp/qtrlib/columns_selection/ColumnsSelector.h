#pragma once

#include <functional>
#include <string>
#include <vector>
#include <filesystem>
#include <future>

#include "FingerprintTable.h"

#include "RawBucketsIO.h"
#include "ColumnsIO.h"

namespace qtr {

    using selection_result_t = std::vector<size_t>;
    using selection_argument_t = const IndigoFingerprintTable &;

    template<typename Functor>
    class ColumnsSelector {
    public:
        ColumnsSelector(std::vector<std::filesystem::path> dataDirectories, Functor selectFunction) :
                _dataDirectories(std::move(dataDirectories)),
                _selectFunction(std::move(selectFunction)) {};

        void handleRawBuckets();

    private:
        std::vector<std::filesystem::path> _dataDirectories;
        Functor _selectFunction;
    };

    static IndigoFingerprintTable readRawBucket(const std::filesystem::path &rawBucketDirPath) {
        IndigoFingerprintTable bucket;
        for (const auto &rawBucketFilePath: findFiles(rawBucketDirPath, rawBucketExtension)) {
            for (const auto &[_, fp]: RawBucketReader(rawBucketFilePath)) {
                bucket.emplace_back(fp);
            }
        }
        return bucket;
    }

    static void saveColumns(const std::vector<size_t> &columns, const std::filesystem::path &rawBucketPath) {
        std::filesystem::path columnsPath = rawBucketPath;
        columnsPath.replace_extension(columnsExtension);
        ColumnsWriter(columnsPath).write(columns);
    }

    template<typename Functor>
    static void handleRawBucket(const std::filesystem::path &rawBucketPath, const Functor &selectFunction) {
        LOG(INFO) << "Start selection columns in " << rawBucketPath;
        auto rawBucket = readRawBucket(rawBucketPath);
        auto selectedColumns = selectFunction(rawBucket);
        saveColumns(selectedColumns, rawBucketPath);
        LOG(INFO) << "Finish selecting columns in " << rawBucketPath;
    }

    template<typename Functor>
    static void handleOneDir(const std::filesystem::path &dirPath, const Functor &selectFunction) {
        auto bucketPaths = findFiles(dirPath, "");
        LOG(INFO) << "Start handling dir: " << dirPath << "with " << bucketPaths.size() << " buckets inside";
        static const size_t step = 32;
        for (size_t i = 0; i < bucketPaths.size(); i += step) {
            std::vector<std::future<void>> tasks;
            for (size_t j = i; j < bucketPaths.size() && j < i + step; j++) {
                tasks.emplace_back(
                        std::async(std::launch::async, handleRawBucket<Functor>, bucketPaths[j], selectFunction)
                );
            }
            for (auto& task : tasks) {
                task.get();
            }
        }
        for (auto &bucketPath: bucketPaths) {
            handleRawBucket(bucketPath, selectFunction);
        }
        LOG(INFO) << "Finish handling dir: " << dirPath << "with " << bucketPaths.size() << " buckets inside";
    }

    template<typename Functor>
    void ColumnsSelector<Functor>::handleRawBuckets() {
        std::vector<std::future<void>> tasks;
        for (auto &dirPath: _dataDirectories) {
            tasks.emplace_back(
                    std::async(std::launch::async, handleOneDir<Functor>, dirPath, _selectFunction)
            );
        }
        for (auto &task: tasks) {
            task.get();
        }
    }

} // namespace qtr
