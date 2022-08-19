#include "SplitterTree.h"
#include "SplitterTreeUtils.h"
#include "Utils.h"

#include <future>
#include <utility>
#include <random>

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

    SplitterTree::Node::Node(SplitterTree *tree, uint64_t depth, std::vector<std::filesystem::path> relatedFiles)
            : _tree(tree), _depth(depth),
              _splitBit(-1),
              _leftChild(nullptr),
              _rightChild(nullptr),
              _relatedFiles(std::move(relatedFiles)) {
        _id = _tree->_countOfNodes++;
    }

    void SplitterTree::build(uint64_t maxDepth, uint64_t maxBucketSize, uint64_t parallelize_depth) const {
        assert(parallelize_depth < maxDepth);
        LOG(INFO) << "Start creating first parallelize_depth(" << parallelize_depth << ") levels of tree";
        auto nodesToRunTasks = _root->buildSubTree(parallelize_depth, maxBucketSize, true);
        std::vector<std::future<std::vector<Node *>>> tasks;
        tasks.reserve(nodesToRunTasks.size());
        LOG(INFO) << "Start sub trees parallelization";
        for (size_t i = 0; i < nodesToRunTasks.size(); i++) {
            size_t baseDirId = std::min(_directories.size() - 1, i / _directories.size());
            tasks.emplace_back(
                    std::async(std::launch::async, &SplitterTree::Node::buildSubTree, nodesToRunTasks[i],
                               maxDepth, maxBucketSize, false, _directories[baseDirId])
            );
        }
        for (auto &task: tasks) {
            task.get();
        }
    }

    void SplitterTree::dump(std::ostream &out) const {
        uint64_t treeSize = size();
        out.write((char *) &treeSize, sizeof treeSize);
        _root->dumpSubTree(out);
    }

    SplitterTree::SplitterTree(std::vector<std::filesystem::path> directories) : _directories(std::move(directories)),
                                                                                 _countOfNodes(0) {
        std::vector<std::filesystem::path> rootRelatedFiles;
        for (auto &dir: _directories) {
            auto rootRelatedDir = dir / "0";
            assert(std::filesystem::is_directory(rootRelatedDir));
            for (auto &filePath: findFiles(rootRelatedDir, ".rb")) {
                rootRelatedFiles.emplace_back(filePath);
            }
        }
        _root = new SplitterTree::Node(this, 0, rootRelatedFiles);
        std::cerr << "root related files: \n";
        for (auto &f: rootRelatedFiles) {
            std::cerr << f << '\n';
        }
    }

    uint64_t SplitterTree::size() const {
        return _countOfNodes;
    }

    std::vector<SplitterTree::Node *>
    SplitterTree::Node::buildSubTree(uint64_t maxDepth, uint64_t maxSizeOfBucket, bool parallelize,
                                     const std::filesystem::path &baseDir) {
        if (_depth >= maxDepth)
            return {this};

        _splitBit = findBestBitToSplit(getFilesPaths(), parallelize);
        auto [leftSize, rightSize] = splitNode(parallelize, baseDir);
        auto leftLeafs = leftSize <= maxSizeOfBucket ? std::vector{_leftChild} :
                         _leftChild->buildSubTree(maxDepth, maxSizeOfBucket, parallelize, baseDir);
        auto rightLeafs = rightSize <= maxSizeOfBucket ? std::vector{_rightChild} :
                          _rightChild->buildSubTree(maxDepth, maxSizeOfBucket, parallelize, baseDir);
        auto leafs = std::move(leftLeafs);
        leafs.insert(leafs.end(), std::move_iterator(rightLeafs.begin()), std::move_iterator(rightLeafs.end()));
        return leafs;
    }

    const std::vector<std::filesystem::path> &SplitterTree::Node::getFilesPaths() const {
        return _relatedFiles;
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

    void SplitterTree::Node::addBucketFilesByParent(SplitterTree::Node *parent) {
        assert(parent != nullptr);
        auto &parentFiles = parent->getFilesPaths();
        for (auto &parentFile: parentFiles) {
            std::filesystem::path baseDir = parentFile.parent_path().parent_path();
            std::filesystem::path currDir = baseDir / std::to_string(_id);
            std::filesystem::create_directory(currDir);
            std::filesystem::path currFile = currDir / parentFile.filename();
            _relatedFiles.emplace_back(currFile);
        }
    }

    void SplitterTree::Node::addBucketFile(const std::filesystem::path &baseDirectory) {
        std::filesystem::path currDir = baseDirectory / std::to_string(_id);
        assert(!std::filesystem::exists(currDir));
        std::filesystem::create_directory(currDir);
        std::filesystem::path currFile = currDir / "data.rb";
        _relatedFiles.emplace_back(currFile);
    }

    void SplitterTree::Node::clearData() {
        for (auto &filePath: _relatedFiles) {
            if (std::filesystem::exists(filePath)) {
                std::filesystem::remove_all(filePath.parent_path());
            }
        }
        _relatedFiles.clear();
    }

    std::pair<uint64_t, uint64_t>
    SplitterTree::Node::splitNode(bool parallelize, const std::filesystem::path &baseDir) {
        assert(_leftChild == nullptr && _rightChild == nullptr && "Node already have been split");
        assert(_splitBit != -1 && "split bit was not defined");

        _leftChild = new SplitterTree::Node(_tree, _depth + 1);
        _rightChild = new SplitterTree::Node(_tree, _depth + 1);
        uint64_t leftSize, rightSize;
        if (parallelize) {
            size_t bucketFilesCount = getFilesPaths().size();
            assert(bucketFilesCount != 0 && "Can not split bucket without files");
            LOG(INFO) << "Start parallel splitting of node " << _id << " with " << bucketFilesCount
                      << " files related to it";
            _leftChild->addBucketFilesByParent(this);
            _rightChild->addBucketFilesByParent(this);
            std::tie(leftSize, rightSize) = splitRawBucketByBitParallel(getFilesPaths(), _splitBit,
                                                                        _leftChild->getFilesPaths(),
                                                                        _rightChild->getFilesPaths());
        } else {
            _leftChild->addBucketFile(baseDir);
            _rightChild->addBucketFile(baseDir);
            assert(_leftChild->getFilesPaths().size() == 1 && _rightChild->getFilesPaths().size() == 1);
            std::tie(leftSize, rightSize) = splitRawBucketByBitNotParallel(getFilesPaths(), _splitBit,
                                                                           _leftChild->getFilesPaths()[0],
                                                                           _rightChild->getFilesPaths()[0]);
        }
        clearData();
        return {leftSize, rightSize};
    }


} // namespace qtr

