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

    using selection_result_t = std::vector<int>;
    using select_argument_t = const IndigoFingerprintTable &;

    template<typename Functor>
    class ColumnsSelector {
    public:
        ColumnsSelector(std::filesystem::path bucketsDir, Functor selectFunction) :
                _bucketsDir(std::move(bucketsDir)),
                _selectFunction(std::move(selectFunction)) {
        };

        void handleRawBuckets();

    private:
        std::filesystem::path _bucketsDir;
        Functor _selectFunction;
    };

    static IndigoFingerprintTable readRawBucket(const std::filesystem::path &rawBucketDirPath) {
        IndigoFingerprintTable bucket;
        for (const auto& rawBucketFilePath : findFiles(rawBucketDirPath, rawBucketExtension)) {
            for (const auto &[_, fp]: RawBucketReader(rawBucketFilePath)) {
                bucket.emplace_back(fp);
            }
        }
        return bucket;
    }

    static void saveColumns(const std::vector<column_t> &columns, const std::filesystem::path &rawBucketPath) {
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
    void ColumnsSelector<Functor>::handleRawBuckets() {
        auto bucketPaths = findFiles(_bucketsDir, "");
        LOG(INFO) << "buckets dir: " << _bucketsDir << " count: " << bucketPaths.size();
        static const size_t step = 20;
        for (size_t i = 0; i < bucketPaths.size(); i += step) {
            std::vector<std::future<void>> tasks;
            for (size_t j = i; j < i + step && j < bucketPaths.size(); j++) {
                const auto& bucketPath = bucketPaths[j];
                tasks.emplace_back(std::async(std::launch::async, handleRawBucket<Functor>, bucketPath, _selectFunction));
            }
            for (auto &task: tasks) {
                task.get();
            }
        }
    }

} // namespace qtr
