#include "SplitterTreeUtils.h"
#include <future>

namespace qtr {
    namespace {
        /**
         * @return {sums for each column in file, number of molecules in file}
        **/
        std::pair<std::vector<uint64_t>, uint64_t>
        findColumnsSumInFile(const std::filesystem::path &rawBucketFilePath) {
            std::vector<uint64_t> columnsSum(IndigoFingerprint::size(), 0);
            uint64_t bucketSize = 0;
            for (const auto &[_, fingerprint]: RawBucketReader(rawBucketFilePath)) {
                bucketSize++;
                for (size_t i = 0; i < IndigoFingerprint::size(); ++i)
                    columnsSum[i] += fingerprint[i];
            }
            return {columnsSum, bucketSize};
        }

        void mergeResults(std::vector<uint64_t> &currColumnsSum, uint64_t &currBucketSize,
                          const std::vector<uint64_t> &fileColumnsSum, uint64_t fileBucketSize) {
            currBucketSize += fileBucketSize;
            for (size_t i = 0; i < IndigoFingerprint::size(); i++) {
                currColumnsSum[i] += fileColumnsSum[i];
            }
        }

        std::pair<std::vector<uint64_t>, uint64_t>
        findColumnsSumInBucketNotParallel(const std::vector<std::filesystem::path> &bucketFiles) {
            std::vector<uint64_t> columnsSum(IndigoFingerprint::size(), 0);
            uint64_t bucketSize = 0;
            for (const auto &filePath: bucketFiles) {
                const auto &[fileColumnsSum, fileBucketSize] = findColumnsSumInFile(filePath);
                mergeResults(columnsSum, bucketSize, fileColumnsSum, fileBucketSize);
            }
            return {columnsSum, bucketSize};
        }

        std::pair<std::vector<uint64_t>, uint64_t>
        findColumnsSumInBucketParallel(const std::vector<std::filesystem::path> &bucketFiles) {
            std::vector<std::future<std::pair<std::vector<uint64_t>, uint64_t>>> tasks;
            tasks.reserve(bucketFiles.size());
            for (const auto &filePath: bucketFiles) {
                tasks.emplace_back(std::async(std::launch::async, findColumnsSumInFile, filePath));
            }
            std::vector<uint64_t> columnsSum(IndigoFingerprint::size(), 0);
            uint64_t bucketSize = 0;
            for (auto &task: tasks) {
                const auto &[fileColumnsSum, fileBucketSize] = task.get();
                mergeResults(columnsSum, bucketSize, fileColumnsSum, fileBucketSize);
            }
            return {columnsSum, bucketSize};
        }

        std::pair<uint64_t, uint64_t>
        splitRawBucketStreamByBit(RawBucketReader &reader, uint64_t &splitBit, RawBucketWriter &zerosWriter,
                                  RawBucketWriter &onesWriter) {
            uint64_t leftSize = 0, rightSize = 0;
            for (const auto &[smiles, fp]: reader) {
                if (fp[splitBit]) {
                    onesWriter.write(make_pair(smiles, fp));
                    leftSize++;
                } else {
                    zerosWriter.write(make_pair(smiles, fp));
                    rightSize++;
                }
            }
            return {leftSize, rightSize};
        }

        std::pair<uint64_t, uint64_t>
        splitFileByBit(const std::filesystem::path &bucketFile, uint64_t splitBit,
                       const std::filesystem::path &zerosFile,
                       const std::filesystem::path &onesFile) {
            RawBucketReader reader(bucketFile);
            RawBucketWriter zerosWriter(zerosFile), onesWriter(onesFile);
            auto [leftSize, rightSize] = splitRawBucketStreamByBit(reader, splitBit, zerosWriter, onesWriter);
            return {leftSize, rightSize};
        }
    } // namespace

    std::pair<uint64_t, uint64_t>
    splitRawBucketByBitNotParallel(const std::vector<std::filesystem::path> &bucketFiles, uint64_t splitBit,
                                   const std::filesystem::path &zeroFilePath,
                                   const std::filesystem::path &onesFilePath) {
        uint64_t leftSize = 0, rightSize = 0;
        RawBucketWriter zerosWriter(zeroFilePath);
        RawBucketWriter onesWriter(onesFilePath);
        for (const auto &filePath: bucketFiles) {
            RawBucketReader reader(filePath);
            auto [fileLeftSize, fileRightSize] = splitRawBucketStreamByBit(reader, splitBit, zerosWriter, onesWriter);
            leftSize += fileLeftSize;
            rightSize += fileRightSize;
        }
        return {leftSize, rightSize};
    }

    std::pair<uint64_t, uint64_t>
    splitRawBucketByBitParallel(const std::vector<std::filesystem::path> &bucketFiles, uint64_t splitBit,
                                const std::vector<std::filesystem::path> &zerosFiles,
                                const std::vector<std::filesystem::path> &onesFiles) {
        assert(bucketFiles.size() == zerosFiles.size() && bucketFiles.size() == onesFiles.size());
        LOG(INFO) << "Start parallel splitting of " << bucketFiles[0] << " and other " << bucketFiles.size() - 1
                  << " neighbor files";
        std::vector<std::future<std::pair<uint64_t, uint64_t>>> tasks;
        tasks.reserve(bucketFiles.size());
        for (size_t i = 0; i < bucketFiles.size(); i++) {
            tasks.emplace_back(std::async(std::launch::async, splitFileByBit, bucketFiles[i], splitBit, zerosFiles[i],
                                          onesFiles[i]));
        }
        uint64_t leftSize = 0, rightSize = 0;
        for (auto &task: tasks) {
            auto [fileLeftSize, fileRightSize] = task.get();
            leftSize += fileLeftSize;
            rightSize += fileRightSize;
        }
        LOG(INFO) << "Finish parallel splitting of " << bucketFiles[0] << " and other " << bucketFiles.size() - 1
                  << " neighbor files";
        return {leftSize, rightSize};
    }

    uint64_t findBestBitToSplit(const std::vector<std::filesystem::path> &rawBucketFiles, bool parallelize) {
        uint64_t bucketSize;
        std::vector<uint64_t> columnsSum;
        if (parallelize) {
            std::tie(columnsSum, bucketSize) = findColumnsSumInBucketParallel(rawBucketFiles);
        } else {
            std::tie(columnsSum, bucketSize) = findColumnsSumInBucketNotParallel(rawBucketFiles);
        }
        uint64_t bestBit = std::min_element(columnsSum.begin(), columnsSum.end(), [bucketSize](uint64_t a, uint64_t b) {
            return std::abs(2 * (int64_t) a - (int64_t) bucketSize) < std::abs(2 * (int64_t) b - (int64_t) bucketSize);
        }) - columnsSum.begin();
        return bestBit;
    }

} // namespace qtr