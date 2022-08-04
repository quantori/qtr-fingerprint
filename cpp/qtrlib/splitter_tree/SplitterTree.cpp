#include "SplitterTree.h"
#include "SplitterTreeUtils.h"

#include <future>

using namespace qtr;

namespace qtr {
    SplitterTree::~SplitterTree() {
        delete _root;
    }

    SplitterTree::Node::~Node() {
        delete _leftChild;
        delete _rightChild;
        std::vector<int> a;
    }

    SplitterTree::Node::Node(SplitterTree *tree, size_t depth) : _tree(tree), _depth(depth), _splitBit(-1),
                                                                 _leftChild(nullptr), _rightChild(nullptr) {
        _id = _tree->_countOfNodes;
    }


    std::vector<std::filesystem::path> SplitterTree::build(size_t maxDepth, uint64_t maxBucketSize) const {
        auto nodesToRunTasks = _root->buildSubTree(PARALLELIZE_DEBT, maxBucketSize);
        using future_t = decltype(std::async(std::launch::async, &SplitterTree::Node::buildSubTree, _root, maxDepth,
                                             maxBucketSize));
        std::vector<future_t> tasks;
        tasks.reserve(nodesToRunTasks.size());
        for (auto leaf: nodesToRunTasks) {
            tasks.emplace_back(
                    std::async(std::launch::async, &SplitterTree::Node::buildSubTree, leaf, maxDepth - PARALLELIZE_DEBT,
                               maxBucketSize)
            );
        }
        std::vector<std::filesystem::path> bucketPaths;
        for (auto &task: tasks) {
            auto leafs = task.get();
            for (auto &leaf: leafs)
                bucketPaths.emplace_back(leaf->getFilePath());
        }
        return bucketPaths;
    }

    void SplitterTree::dump(std::ostream &out) const {
        size_t treeSize = size();
        out.write((char *) &treeSize, sizeof treeSize);
        _root->dumpSubTree(out);
    }

    SplitterTree::SplitterTree(std::filesystem::path directory) : _directory(std::move(directory)),
                                                                  _root(new Node(this, 0)), _countOfNodes(0) {}

    std::vector<SplitterTree::Node *>
    SplitterTree::Node::buildSubTree(size_t maxDepth, uint64_t maxSizeOfBucket) {
        if (_depth >= maxDepth)
            return {this};

        _splitBit = findBestBitToSplit(getFilePath());
        _leftChild = new SplitterTree::Node(_tree, _depth + 1);
        _rightChild = new SplitterTree::Node(_tree, _depth + 1);
        auto [leftSize, rightSize] = splitRawBucketByBit(getFilePath(), _splitBit, _leftChild->getFilePath(),
                                                         _rightChild->getFilePath());
        auto leftLeafs = leftSize <= maxSizeOfBucket ? std::vector{_leftChild} :
                         _leftChild->buildSubTree(maxDepth, maxSizeOfBucket);
        auto rightLeafs = rightSize <= maxSizeOfBucket ? std::vector{_rightChild} :
                          _rightChild->buildSubTree(maxDepth, maxSizeOfBucket);
        auto leafs = std::move(leftLeafs);
        leafs.insert(leafs.end(), std::move_iterator(rightLeafs.begin()), std::move_iterator(rightLeafs.end()));
        return leafs;
    }

    std::filesystem::path SplitterTree::Node::getFilePath() const {
        return _tree->_directory / std::to_string(_id);
    }

    void SplitterTree::Node::dumpSubTree(std::ostream &out) {
        out.write((char *) &_id, sizeof _id);
        out.write((char *) &_splitBit, sizeof _splitBit);
        uint64_t leftId = SplitterTree::Node::getId(_leftChild);
        uint64_t rightId = SplitterTree::Node::getId(_rightChild);
        out.write((char *) &leftId, sizeof(leftId));
        out.write((char *) &rightId, sizeof(rightId));

        if (_leftChild != nullptr)
            _leftChild->dumpSubTree(out);
        if (_rightChild != nullptr)
            _rightChild->dumpSubTree(out);
    }

    uint64_t SplitterTree::Node::getId(SplitterTree::Node *node) {
        return node == nullptr ? -1 : node->_id;
    }

} // namespace qtr

