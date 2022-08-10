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

    using choice_result_t = std::vector<int>;
    using choice_argumnent_t = const IndigoFingerprintTable &;

    template<typename Functor>
    class ColumnsChooser {
    public:
        ColumnsChooser(std::filesystem::path bucketsDir, Functor chooseFunc) :
                _bucketsDir(std::move(bucketsDir)),
                _choiceFunc(std::move(chooseFunc)) {
        };

        void handleRawBuckets();

    private:
        std::filesystem::path _bucketsDir;
        Functor _choiceFunc;
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
    static void handleRawBucket(const std::filesystem::path &rawBucketPath, const Functor &choiceFunc) {
        auto rawBucket = readRawBucket(rawBucketPath);
        auto chosenColumns = choiceFunc(rawBucket);
        saveColumns(chosenColumns, rawBucketPath);
    }

    template<typename Functor>
    void ColumnsChooser<Functor>::handleRawBuckets() {
        using future_t = decltype(std::async(std::launch::async, handleRawBucket<Functor>, "", _choiceFunc));
        std::vector<future_t> tasks;
        int started = 0;
        for (const auto &bucketPath: findFiles(_bucketsDir, "")) {
            assert(std::filesystem::is_directory(bucketPath));
            tasks.emplace_back(std::async(std::launch::async, handleRawBucket<Functor>, bucketPath, _choiceFunc));
            LOG(INFO) << "Start choosing columns for " << bucketPath << " (" << ++started << ")";
        }
        int completed = 0;
        for (auto &task: tasks) {
            task.get();
            LOG(INFO) << "Complete choosing columns (" << ++completed << ")";
        }
    }

} // namespace qtr
