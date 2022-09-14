#pragma once

#include <utility>
#include <cstdint>
#include <filesystem>

#include "SplitterTree.h"
#include "raw_bucket_io/RawBucketReader.h"
#include "raw_bucket_io/RawBucketWriter.h"

namespace qtr {
    /**
     * Finds column, that is the best to split with.
     * Criteria: absolute difference between records with ones in that column and records with zeroes, is minimum.
     * @return column id
     */
    uint64_t findBestBitToSplit(const std::vector<std::filesystem::path> &rawBucketFiles, bool parallelize);


    std::pair<uint64_t, uint64_t>
    splitRawBucketByBitParallel(const std::vector<std::filesystem::path> &bucketFiles, uint64_t splitBit,
                                const std::vector<std::filesystem::path> &zerosFiles,
                                const std::vector<std::filesystem::path> &onesFiles);

    std::pair<uint64_t, uint64_t>
    splitRawBucketByBitNotParallel(const std::vector<std::filesystem::path> &bucketFiles, uint64_t splitBit,
                                   const std::filesystem::path &zeroFilePath,
                                   const std::filesystem::path &onesFilePath);
}