#pragma once

#include <utility>
#include <cstdint>
#include <filesystem>

#include "SplitterTree.h"
#include "RawBucketsIO.h"

namespace qtr {

    /**
     * @brief Split given bucket in two buckets by value in given column
     * @param fileToSplitPath path to file which must be split
     * @b Note: this file is deleted after function execution
     * @param splitBit index of the column by the values in which bucket must be split
     * @param zeroDirPath path to file where to store bucket with zeros on @c splitBit column
     * @param onesDirPath path to file where to store bucket with ones on @c splitBit column
     * @return <count of records in bucket with zeros, count of records in bucket with ones>
     */
    std::pair<uint64_t, uint64_t> splitRawBucketByBit(const std::filesystem::path &fileToSplitPath, uint64_t splitBit,
                                                      const std::filesystem::path &zeroDirPath,
                                                      const std::filesystem::path &onesDirPath,
                                                      bool parallelizeByFiles);


    /**
     * Finds column, that is the best to split with.
     * Criteria: absolute difference between records with ones in that column and records with zeroes, is minimum.
     * @return column id
     */
    uint64_t findBestBitToSplit(const std::filesystem::path &rawBucketPath);


    std::vector<std::filesystem::path> nodesToFilePaths(const std::vector<SplitterTree::Node *> &nodes);
}