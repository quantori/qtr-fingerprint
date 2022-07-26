#pragma once

#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <cstdio>
#include <future>

#include "Fingerprint.h"
#include "Utils.h"

namespace qtr {
    const int COUNT_THREADS_POW = 4; // should be equals to log2(your count of threads)

    class SplitterTree {
        static uint64_t _countOfNodes; // default set to zero, increases after each constructor call
        uint64_t _nodeNumber;
        std::string _filename;
        size_t _depth;

        uint64_t _splitBit;
        SplitterTree *_leftChild = nullptr; // go left if bit is zero
        SplitterTree *_rightChild = nullptr; // go right if bit is one
    public:
        /**
         * Not building tree, just save data for future splitting.
         * @param filename name of file with data
         * @param depth depth of node
         */
        inline SplitterTree(std::string filename, size_t depth = 0) : _filename(filename), _depth(depth) {
            _nodeNumber = _countOfNodes;
            _countOfNodes++;
        }

        /**
         * Building a tree with given depth. This node should be a leaf before call. Otherwise UB
         * After splitting, current file(_filename) will be deleted.
         * @param maxDepth max depth of tree
         * @param minCountInNode if count of records is less than [minCountInNode], than method do nothing
         * @return vector of filenames(buckets), size of vector is not more than 2^maxDepth.
         */
        std::vector<std::string> split(size_t maxDepth, uint64_t minCountInNode);

        /**
         * Saves tree with root at this node to stream
         * Format of saving: _nodeNumber|_splitBit|_leftChild->number|_rightChild->number
         * Like 4 unsigned integral types in a row, without "|", if _leftChild is null,
         * than we take it's number as (uint64_t)-1
         * @param out
         */
        void saveTo(std::ostream &out);

    private:
        /**
         * @param node
         * @return -1 if node is nullptr and node number otherwise
         */
        uint64_t getNum(SplitterTree *node);
    };

}
