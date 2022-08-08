#include "SplitterTree.h"
#include "SplitterTreeUtils.h"
#include "Utils.h"

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

    SplitterTree::Node::Node(SplitterTree *tree, uint64_t depth) : _tree(tree), _depth(depth),
                                                                   _splitBit(-1),
                                                                   _leftChild(nullptr),
                                                                   _rightChild(nullptr),
                                                                   _filesNumber(0) {
        _id = _tree->_countOfNodes++;
        std::filesystem::create_directory(getDirPath());
        _filesNumber = findFiles(getDirPath(), rawBucketExtension).size();
    }

    std::vector<std::filesystem::path>
    SplitterTree::buildWithoutSubTreeParallelize(uint64_t maxDepth, uint64_t maxBucketSize) const {
        auto leafs = _root->buildSubTree(maxDepth, maxBucketSize, true);
        return nodesToDirPaths(leafs);
    }

    std::vector<std::filesystem::path>
    SplitterTree::buildWithSubTreeParallelize(uint64_t maxDepth, uint64_t maxBucketSize,
                                              uint64_t parallelize_depth) const {
        assert(parallelize_depth < maxDepth);
        auto nodesToRunTasks = _root->buildSubTree(parallelize_depth, maxBucketSize, true);
        std::vector<std::future<std::vector<Node *>>> tasks;
        tasks.reserve(nodesToRunTasks.size());
        for (auto leaf: nodesToRunTasks) {
            tasks.emplace_back(
                    std::async(std::launch::async, &SplitterTree::Node::buildSubTree, leaf,
                               maxDepth - parallelize_depth, maxBucketSize, false)
            );
        }
        std::vector<std::filesystem::path> bucketPaths;
        for (auto &task: tasks) {
            auto leafs = task.get();
            auto directoriesPaths = nodesToDirPaths(leafs);
            bucketPaths.insert(bucketPaths.end(), std::make_move_iterator(directoriesPaths.begin()),
                               std::make_move_iterator(directoriesPaths.end()));
        }
        return bucketPaths;
    }

    std::vector<std::filesystem::path>
    SplitterTree::build(uint64_t maxDepth, uint64_t maxBucketSize, uint64_t parallelize_depth) const {
        if (maxDepth <= parallelize_depth) {
            return buildWithoutSubTreeParallelize(maxDepth, maxBucketSize);
        } else {
            return buildWithSubTreeParallelize(maxDepth, maxBucketSize, parallelize_depth);
        }
    }

    void SplitterTree::dump(std::ostream &out) const {
        uint64_t treeSize = size();
        out.write((char *) &treeSize, sizeof treeSize);
        _root->dumpSubTree(out);
    }

    SplitterTree::SplitterTree(std::filesystem::path directory) : _directory(std::move(directory)), _countOfNodes(0),
                                                                  _root(new Node(this, 0)) {}

    uint64_t SplitterTree::size() const {
        return _countOfNodes;
    }

    std::vector<SplitterTree::Node *>
    SplitterTree::Node::buildSubTree(uint64_t maxDepth, uint64_t maxSizeOfBucket, bool parallelize) {
        if (_depth >= maxDepth)
            return {this};

        _splitBit = findBestBitToSplit(getDirPath(), parallelize);
        _leftChild = new SplitterTree::Node(_tree, _depth + 1);
        _rightChild = new SplitterTree::Node(_tree, _depth + 1);
        uint64_t leftSize, rightSize;
        if (parallelize) {
            size_t bucketFilesCount = getFilesPaths().size();
            _leftChild->addBucketFiles(bucketFilesCount);
            _rightChild->addBucketFiles(bucketFilesCount);
            std::tie(leftSize, rightSize) = splitRawBucketByBitParallel(getFilesPaths(), _splitBit,
                                                                        _leftChild->getFilesPaths(),
                                                                        _rightChild->getFilesPaths());
        } else {
            _leftChild->addBucketFiles(1);
            _rightChild->addBucketFiles(1);
            std::tie(leftSize, rightSize) = splitRawBucketByBitNotParallel(getFilesPaths(), _splitBit,
                                                                           _leftChild->getFilesPaths()[0],
                                                                           _rightChild->getFilesPaths()[0]);
        }
        auto leftLeafs = leftSize <= maxSizeOfBucket ? std::vector{_leftChild} :
                         _leftChild->buildSubTree(maxDepth, maxSizeOfBucket, parallelize);
        auto rightLeafs = rightSize <= maxSizeOfBucket ? std::vector{_rightChild} :
                          _rightChild->buildSubTree(maxDepth, maxSizeOfBucket, parallelize);
        auto leafs = std::move(leftLeafs);
        leafs.insert(leafs.end(), std::move_iterator(rightLeafs.begin()), std::move_iterator(rightLeafs.end()));
        return leafs;
    }

    // todo refactor of places of usage
    std::filesystem::path SplitterTree::Node::getDirPath() const {
        return _tree->_directory / std::to_string(_id);
    }

    std::vector<std::filesystem::path> SplitterTree::Node::getFilesPaths() const {
        auto paths = findFiles(getDirPath(), rawBucketExtension);
        assert(paths.size() == _filesNumber);
        return paths;
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

    void SplitterTree::Node::addBucketFiles(uint64_t count) {
        auto dirPath = getDirPath();
        for (size_t i = 0; i < count; i++) {
            auto filePath = dirPath / (std::to_string(_filesNumber));
            filePath.replace_extension(rawBucketExtension);
            assert(!std::filesystem::exists(filePath) && "bucket already exists");
            std::ofstream fileCreator(filePath);
            _filesNumber++;
        }
    }

} // namespace qtr

