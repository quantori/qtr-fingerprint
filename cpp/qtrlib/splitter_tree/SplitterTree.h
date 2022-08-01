#pragma once

#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <cstdio>
#include <future>
#include <algorithm>
#include <cmath>

#include "Fingerprint.h"
#include "Utils.h"

namespace qtr {
    const int COUNT_THREADS_POW = 3; // should be equals to log2(your count of threads)

    /**
     * for better understanding formats of files, see
     * https://quantori.atlassian.net/wiki/spaces/QLAB/pages/3139567644/1st+implementation+ideas#%D0%A4%D0%BE%D1%80%D0%BC%D0%B0%D1%82%D1%8B-%D1%84%D0%B0%D0%B9%D0%BB%D0%BE%D0%B2%3A
     */
    class SplitterTree {
        static std::atomic_uint64_t _countOfNodes; // default set to zero, increases after each constructor call
        uint64_t _nodeNumber;
        std::string _filename;
        std::filesystem::path _dir; // location dir
        size_t _depth;

        uint64_t _splitBit = -1;
        SplitterTree *_leftChild = nullptr; // go left if bit is zero
        SplitterTree *_rightChild = nullptr; // go right if bit is one
    public:
        inline ~SplitterTree() {
            delete _leftChild;
            delete _rightChild;
        }

        /**
         * Not building tree, just save data for future splitting.
         * @param dir directory, so all file names will be "dir / filename"
         * @param filename name of file with data
         * @param depth depth of node
         */
        inline SplitterTree(const std::filesystem::path &dir, const std::string &filename, size_t depth = 0) :
                _dir(dir), _filename(filename), _depth(depth) {
            _nodeNumber = _countOfNodes++;
        }

        /**
         * Building a tree with given depth. This node should be a leaf before call. Otherwise UB
         * After splitting, current file(_filename) will be deleted
         * @param maxDepth max depth of tree
         * @param minCountInNode if count of records is less than [minCountInNode], than method do nothing
         * @return vector of filenames(buckets), size of vector is not more than 2^maxDepth.
         */
        std::vector<std::string> split(size_t maxDepth, uint64_t minCountInNode);

        /**
         * Saves tree with root at this node to stream
         * Format of saving:
         * count of nodes in tree, and then each node like 4 unsigned integral types in a row, without "|",
         * _nodeNumber|_splitBit|_leftChild->number|_rightChild->number
         *
         * if _leftChild is null,x`
         * than we take it's number as (uint64_t)-1
         * @param out
         * @param writeSize if false, than count of nodes, wont be wrote
         */
        void saveTo(std::ostream &out, bool writeSize = true);

        // TODO delete unused functions connected with correlated columns

        /**
         * Saves mostly non correlating order of columns for each leaf, file format:
         * filename: leaf->_filenameOrderColumns, for example: "1OrderColumns"
         * and inside each file we have an order of columns, size_t numbers split with a space
         * Example: "5 2 3 1 4 0"
         */
        void saveNonCorrelatingColumnInEachBucket();

        /**
         * @return count of nodes in a tree
         */
        uint64_t size() const;

    private:
        /**
         * @param node
         * @return -1 if node is nullptr and node number otherwise
         */
        uint64_t getNum(SplitterTree *node);

        /**
         * Finds column, that is the best to split with.
         * Criteria: absolute difference between records with ones in that column and records with zeroes, is minimum.
         * @return column id
         */
        uint64_t findBestBitToSplit();

        /**
         * Deletes current file, and creates two files for children
         * To left child we write records, when bitSplit is zero, and to right child when bitSplit is one
         * @return <count of records in left child, count of records in right child>
         */
        std::pair<uint64_t, uint64_t> prepareFilesForChildren(uint64_t bitSplit);

        /**
         * Runs in children and then merging their results
         * @param sizeLeft - count of records in left child
         * @param sizeRight - count of records in right child
         * @param maxDepth
         * @param minCountInNode
         * @return merged results of children splits
         */
        std::vector<std::string>
        splitInChildren(uint64_t sizeLeft, uint64_t sizeRight, size_t maxDepth, uint64_t minCountInNode);
    };

}
