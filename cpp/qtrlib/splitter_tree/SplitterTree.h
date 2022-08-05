#pragma once

#include <cstring>
#include <string>
#include <utility>
#include <vector>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <cmath>

#include "Fingerprint.h"
#include "Utils.h"

namespace qtr {
    /**
     * @brief <a href="https://quantori.atlassian.net/wiki/spaces/QLAB/pages/3139567644/1st+implementation+ideas#%D0%A4%D0%BE%D1%80%D0%BC%D0%B0%D1%82%D1%8B-%D1%84%D0%B0%D0%B9%D0%BB%D0%BE%D0%B2%3A">There</a>
     * is files format documentation
     *
     */
    class SplitterTree {
    public:
        class Node;

    private:
        std::filesystem::path _directory;
        std::atomic_uint64_t _countOfNodes;
        Node *_root;

        std::vector<std::filesystem::path> buildNotParallel(uint64_t maxDepth, uint64_t maxBucketSize) const;

        std::vector<std::filesystem::path>
        buildParallel(uint64_t maxDepth, uint64_t maxBucketSize, uint64_t parallelize_depth) const;

    public:
        /**
         * @brief Saves data for future building. Doesn't build the tree
         * @param directory directory to save tree's data
         */
        explicit SplitterTree(std::filesystem::path directory);

        ~SplitterTree();

        /**
         * @brief Builds a tree with given parameters.
         * If call this function more than once, you get UB
         * @param maxDepth max depth of tree
         * @param maxBucketSize max number of elements in a leaf
         * @return vector of buckets' filenames
         */
        std::vector<std::filesystem::path>
        build(uint64_t maxDepth, uint64_t maxBucketSize, uint64_t parallelize_debt) const;

        /**
         * @return count of nodes in a tree
         */
        [[nodiscard]] uint64_t size() const;

        /**
         * @brief Dumps tree to @c out.
         *
         * @File_Format
         * 1st 4 bytes are @c uint32_t with count of nodes in tree.
         * Whole next file consist of 16 bytes blocks. Each block consist of 4 @c uint32_t numbers describing node:
         * @c nodeId, @c splitBit, @c leftChildId, @c rightChildId
         *
         * @b Note: if node is leaf than splitBit, leftChildId, rightChildId are equal to -1
         * @param out
         */
        void dump(std::ostream &out) const;
    };

    class SplitterTree::Node {
    public:
        /**
         * @brief Saves data for future building. Doesn't build subtree
         * @param tree tree to which this node belongs
         * @param depth depth of node
         */
        explicit Node(SplitterTree *tree, uint64_t depth);

        ~Node();

        /**
         * @return path to file to save node's data
         */
        [[nodiscard]] std::filesystem::path getFilePath() const;

        /**
         * @param node
         * @return -1 if node is nullptr and node's id otherwise
         */
        static uint64_t getId(Node *node);

        std::vector<SplitterTree::Node *> buildSubTree(uint64_t maxDepth, uint64_t maxSizeOfBucket);

        void dumpSubTree(std::ostream &out);

    private:
        SplitterTree *_tree;
        uint64_t _depth;
        uint64_t _splitBit;
        Node *_leftChild; // go left if bit is zero
        Node *_rightChild; // go right if bit is one
        uint64_t _id;
    };

} // namespace qtr
