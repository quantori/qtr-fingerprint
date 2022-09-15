#include <utility>
#include <future>
#include <cstdlib>
#include <random>

#include "BallTree.h"
#include "data_io/fingerprint_table_io/FingerprintTableReader.h"
#include "data_io/fingerprint_table_io/FingerprintTableWriter.h"

namespace qtr {

    namespace {

        const std::string dataFilename = "data.ft";
        std::mt19937 randomGenerator(0);

        void splitFileBybit(const std::filesystem::path &file, size_t splitBit, const std::filesystem::path &leftFile,
                            const std::filesystem::path &rightFile, size_t leftSize, size_t rightSize) {
            FingerprintTableWriter leftWriter(leftFile), rightWriter(rightFile);
            for (const auto &value: FingerprintTableReader(file)) {
                const auto &[_, fingerprint] = value;
                if (leftSize != 0 && (!fingerprint[splitBit] || rightSize == 0)) {
                    leftWriter << value;
                    leftSize--;
                } else if (rightSize != 0 && (fingerprint[splitBit] || leftSize == 0)) {
                    rightWriter << value;
                    rightSize--;
                } else {
                    assert(false);
                }
            }
            assert(leftSize == 0 && rightSize == 0);
            std::filesystem::remove(file);
        }

        void prepareManyDirsForNode(const std::vector<std::filesystem::path> &dataDirectories, size_t nodeId) {
            for (auto &dirPath: dataDirectories) {
                std::filesystem::create_directory(dirPath / std::to_string(nodeId));
            }
        }

        std::vector<std::filesystem::path>
        findDrivesForNodes(size_t nodesNumber, std::vector<std::filesystem::path> dataDirectories) {
            std::shuffle(dataDirectories.begin(), dataDirectories.end(), randomGenerator);
            std::vector<std::filesystem::path> result(nodesNumber);
            for (size_t i = 0; i < nodesNumber; i++) {
                result[i] = dataDirectories[i % dataDirectories.size()];
            }
            return result;
        }

        ColumnsStatistic createColumnsStatistic(const std::filesystem::path &filePath) {
            return ColumnsStatistic(filePath);
        }

        std::pair<std::vector<ColumnsStatistic>, ColumnsStatistic>
        collectFilesStatistic(const std::vector<std::filesystem::path> &filePaths) {
            std::vector<std::future<ColumnsStatistic>> tasks;
            for (auto &filePath: filePaths) {
                tasks.emplace_back(std::async(std::launch::async, createColumnsStatistic, filePath));
            }
            std::vector<ColumnsStatistic> filesStat;
            ColumnsStatistic summaryStat;
            for (auto &task: tasks) {
                filesStat.emplace_back(task.get());
                summaryStat += filesStat.back();
            }
            return {filesStat, summaryStat};
        }

        std::vector<std::pair<size_t, size_t>>
        findSplitPolicy(const std::vector<ColumnsStatistic> &filesStats, const ColumnsStatistic &summaryStat,
                        size_t splitBit) {
            size_t leftSize = summaryStat.size() / 2;
            size_t rightSize = (summaryStat.size() + 1) / 2;
            std::vector<std::pair<size_t, size_t>> result;
            for (const auto &filesStat: filesStats) {
                size_t fileZeros = filesStat.zeros(splitBit);
                size_t fileOnes = filesStat.ones(splitBit);
                assert(fileZeros + fileOnes <= leftSize + rightSize);
                if (fileZeros <= leftSize && fileOnes <= rightSize) {
                    result.emplace_back(fileZeros, fileOnes);
                } else if (fileZeros > leftSize && fileOnes <= rightSize) {
                    assert(fileOnes + fileZeros >= leftSize);
                    result.emplace_back(leftSize, fileZeros + fileOnes - leftSize);
                } else if (fileZeros <= leftSize && fileOnes > rightSize) {
                    assert(fileZeros + fileOnes >= rightSize);
                    result.emplace_back(fileZeros + fileOnes - rightSize, rightSize);
                } else {
                    assert(false);
                }
                leftSize -= result.back().first;
                rightSize -= result.back().second;
            }
            assert(leftSize == 0 && rightSize == 0);
            return result;
        }

        std::pair<size_t, size_t> findSplitPolicy(const ColumnsStatistic &statistic, size_t splitBit) {
            return {statistic.size() / 2, (statistic.size() + 1) / 2};
        }

        std::filesystem::path switchFilesNode(const std::filesystem::path &filePath, size_t newNodeId) {
            return filePath.parent_path().parent_path() / std::to_string(newNodeId) / filePath.filename();
        }

        void mergeFiles(const std::vector<std::filesystem::path> &sources, const std::filesystem::path &destination) {
            FingerprintTableWriter writer(destination);
            for (auto &source: sources) {
                FingerprintTableReader reader(source);
                std::copy(reader.begin(), reader.end(), writer.begin());
            }
            for (auto &source: sources) {
                std::filesystem::remove(source);
            }
        }

        void putAnswer(std::vector<size_t> &result, size_t value, size_t ansCount, std::mutex &resultLock,
                       bool &isTerminate) {
            std::lock_guard<std::mutex> lock(resultLock);
            result.emplace_back(value);
            isTerminate |= result.size() >= ansCount;
        }

        void searchInFile(const std::filesystem::path &filePath, const IndigoFingerprint &query, size_t ansCount,
                          std::vector<size_t> &result, std::mutex &resultsLock, bool &isTerminate) {
            for (const auto &[id, fingerprint]: FingerprintTableReader(filePath)) {
                if (query <= fingerprint)
                    putAnswer(result, id, ansCount, resultsLock, isTerminate);
            }
        }

    } // namespace

    BallTree::BallTree(size_t depth, size_t parallelizationDepth, std::vector<std::filesystem::path> dataDirectories,
                       const BitSelector &bitSelector)
            : _dataDirectories(std::move(dataDirectories)), _depth(depth) {
        _nodes.resize((1ull << (depth + 1)) - 1);
        _leafDataPaths.resize(1ull << depth);
        buildTree(depth, parallelizationDepth, bitSelector);
    }

    size_t BallTree::leftChild(size_t nodeId) {
        return nodeId * 2 + 1;
    }

    size_t BallTree::rightChild(size_t nodeId) {
        return nodeId * 2 + 2;
    }

    size_t BallTree::root() const {
        return 0;
    }

    size_t BallTree::parent(size_t nodeId) {
        return (nodeId - 1) >> 1ull;
    }

    bool BallTree::isLeaf(size_t nodeId) const {
        return leftChild(nodeId) >= _nodes.size();
    }

    void BallTree::buildTree(size_t depth, size_t parallelizationDepth, const BitSelector &bitSelector) {
        assert(depth > parallelizationDepth);
        LOG(INFO) << "Start tree building";
        LOG(INFO) << "Start first levels building";
        std::vector<size_t> nodes = buildFirstLevels(root(), parallelizationDepth, bitSelector);
        LOG(INFO) << "Finish first levels building";
        std::shuffle(nodes.begin(), nodes.end(), randomGenerator);
        std::vector<std::filesystem::path> nodesDirs = findDrivesForNodes(nodes.size(), _dataDirectories);
        LOG(INFO) << "Start subtree parallelization preparation";
        for (size_t i = 0; i < nodes.size(); i++) {
            std::filesystem::path destinationPath = nodesDirs[i] / std::to_string(nodes[i]) / dataFilename;
            LOG(INFO) << "Start transferring node " << nodes[i] << " to " << destinationPath;
            std::filesystem::create_directory(destinationPath.parent_path());
            std::vector<std::filesystem::path> nodeFiles = getNodeFiles(nodes[i]);
            mergeFiles(nodeFiles, destinationPath);
            LOG(INFO) << "Finish transferring node " << nodes[i] << " to " << destinationPath;
            for (auto &nodeFile: nodeFiles) {
                if (nodeFile.parent_path() != destinationPath.parent_path() &&
                    std::filesystem::exists(nodeFile.parent_path())) {
                    LOG(INFO) << "Delete " << nodeFile.parent_path();
                    std::filesystem::remove_all(nodeFile.parent_path());
                }
            }
        }
        LOG(INFO) << "Finish subtree parallelization preparation";
        LOG(INFO) << "Start last levels building";
        std::vector<std::future<void>> tasks;
        for (size_t i = 0; i < nodes.size(); i++) {
            tasks.emplace_back(std::async(std::launch::async, &BallTree::buildLastLevels, this, nodes[i],
                                          depth - parallelizationDepth, std::cref(bitSelector), nodesDirs[i])
            );
        }
        for (auto &task: tasks) {
            task.get();
        }
        LOG(INFO) << "Finish last levels building";
        LOG(INFO) << "Start centroids calculating";
        calculateCentroid(root());
        LOG(INFO) << "Finish centroids calculating";
        LOG(INFO) << "Finish tree building";
    }

    void BallTree::splitNodeManyFiles(size_t nodeId, const BitSelector &bitSelector) {
        std::vector<std::filesystem::path> nodeFiles = getNodeFiles(nodeId);
        LOG(INFO) << "Start parallel splitting of node " << nodeId << " with " << nodeFiles.size() << " files";
        std::shuffle(nodeFiles.begin(), nodeFiles.end(), randomGenerator);
        LOG(INFO) << "Start collecting statistics for node " << nodeId;
        auto [filesStatistics, summaryStat] = collectFilesStatistic(nodeFiles);
        LOG(INFO) << "Finish collecting statistics for node " << nodeId;
        LOG(INFO) << "Start selecting split bit for node " << nodeId;
        size_t splitBit = bitSelector(summaryStat);
        LOG(INFO) << "Finish selecting split bit for node " << nodeId << ". Selected bit - " << splitBit
                  << ", ones - " << summaryStat.ones(splitBit) << ", zeros - " << summaryStat.zeros(splitBit);

        std::vector<std::pair<size_t, size_t>> splitPolicy = findSplitPolicy(filesStatistics, summaryStat, splitBit);
        prepareManyDirsForNode(_dataDirectories, leftChild(nodeId));
        prepareManyDirsForNode(_dataDirectories, rightChild(nodeId));
        std::vector<std::future<void>> tasks;
        for (size_t i = 0; i < nodeFiles.size(); i++) {
            std::filesystem::path leftFile = switchFilesNode(nodeFiles[i], leftChild(nodeId));
            std::filesystem::path rightFile = switchFilesNode(nodeFiles[i], rightChild(nodeId));
            tasks.emplace_back(
                    std::async(std::launch::async, splitFileBybit, nodeFiles[i], splitBit, leftFile, rightFile,
                               splitPolicy[i].first, splitPolicy[i].second)
            );
        }
        for (auto &task: tasks) {
            task.get();
        }
        LOG(INFO) << "Finish parallel splitting of node " << nodeId;
        for (auto &nodeFile: nodeFiles) {
            if (std::filesystem::exists(nodeFile.parent_path())) {
                LOG(INFO) << "Delete " << nodeFile.parent_path();
                std::filesystem::remove_all(nodeFile.parent_path());
            }
        }
    }

    void BallTree::splitNodeOneFile(size_t nodeId, const BitSelector &bitSelector,
                                    const std::filesystem::path &dataDirectory) {
        std::filesystem::path nodeFile = dataDirectory / std::to_string(nodeId) / dataFilename;
        auto [_, stat] = collectFilesStatistic({nodeFile});
        size_t splitBit = bitSelector(stat);
        auto [leftSize, rightSize] = findSplitPolicy(stat, splitBit);
        std::filesystem::path leftFile = dataDirectory / std::to_string(leftChild(nodeId)) / dataFilename;
        std::filesystem::path rightFile = dataDirectory / std::to_string(rightChild(nodeId)) / dataFilename;
        std::filesystem::create_directory(leftFile.parent_path());
        std::filesystem::create_directory(rightFile.parent_path());
        splitFileBybit(nodeFile, splitBit, leftFile, rightFile, leftSize, rightSize);
        std::filesystem::remove_all(nodeFile.parent_path());
    }


    std::vector<size_t> BallTree::buildFirstLevels(size_t nodeId, size_t levels, const BitSelector &bitSelector) {
        if (levels == 0)
            return {nodeId};
        splitNodeManyFiles(nodeId, bitSelector);
        std::vector<size_t> nodes = buildFirstLevels(leftChild(nodeId), levels - 1, bitSelector);
        std::vector<size_t> rightNodes = buildFirstLevels(rightChild(nodeId), levels - 1, bitSelector);
        nodes.insert(nodes.end(), rightNodes.begin(), rightNodes.end());
        return nodes;
    }

    std::vector<std::filesystem::path> BallTree::getNodeFiles(size_t nodeId) const {
        std::vector<std::filesystem::path> result;
        for (auto &dirPath: _dataDirectories) {
            if (!std::filesystem::exists(dirPath / std::to_string(nodeId)))
                continue;
            for (auto &filePath: findFiles(dirPath / std::to_string(nodeId), ".ft")) {
                result.emplace_back(filePath);
            }
        }
        return result;
    }

    void BallTree::buildLastLevels(size_t nodeId, size_t levels, const BitSelector &bitSelector,
                                   const std::filesystem::path &dataDirectory) {
        if (levels == 0) {
            assert(isLeaf(nodeId));
            return;
        }
        assert(!isLeaf(nodeId));
        splitNodeOneFile(nodeId, bitSelector, dataDirectory);
        buildLastLevels(leftChild(nodeId), levels - 1, bitSelector, dataDirectory);
        buildLastLevels(rightChild(nodeId), levels - 1, bitSelector, dataDirectory);
    }

    void BallTree::calculateCentroid(size_t nodeId) {
        if (isLeaf(nodeId)) {
            _nodes[nodeId].centroid.reset();
            std::vector<std::filesystem::path> nodeFiles = getNodeFiles(nodeId);
            assert(nodeFiles.size() == 1);
            for (const auto &[_, fingerprint]: FingerprintTableReader(nodeFiles[0])) {
                _nodes[nodeId].centroid |= fingerprint;
            }
        } else {
            calculateCentroid(leftChild(nodeId));
            calculateCentroid(rightChild(nodeId));
            _nodes[nodeId].centroid = _nodes[leftChild(nodeId)].centroid;
            _nodes[nodeId].centroid |= _nodes[rightChild(nodeId)].centroid;
        }
    }

    void BallTree::initLeafDataPaths() {
        size_t expectedFilesNumber = (_nodes.size() + 1) / 2;
        _leafDataPaths.resize(expectedFilesNumber);
        std::vector<bool> isInit(expectedFilesNumber, false);
        for (auto &dataDir: _dataDirectories) {
            for (auto &filePath: findFiles(dataDir, "")) {
                size_t nodeId = atoll(filePath.filename().c_str());
                assert(!isInit[nodeId]);
                isInit[nodeId] = true;
            }
        }
        assert(std::count(isInit.begin(), isInit.end(), false) == 0);
    }

    std::vector<size_t> BallTree::search(const IndigoFingerprint &query, size_t ansCount, size_t startDepth) {
        std::vector<std::future<void>> tasks;
        bool isTerminate = false;
        std::vector<size_t> results;
        std::mutex resultsLock;
        for (size_t i = (1ull << startDepth) - 1; i < (1ull << (startDepth - 1)) - 1; i++) {
            tasks.emplace_back(
                    std::async(std::launch::async, &BallTree::searchInSubtree, this, i, std::cref(query), ansCount,
                               std::ref(results), std::ref(resultsLock), std::ref(isTerminate))
            );
        }
        for (auto &task: tasks) {
            task.get();
        }
        return results;
    }

    void BallTree::searchInSubtree(size_t nodeId, const IndigoFingerprint &query, size_t ansCount,
                                   std::vector<size_t> &result,
                                   std::mutex &resultLock, bool &isTerminate) {
        if (isTerminate || !(query <= _nodes[nodeId].centroid)) {
            return;
        }
        if (isLeaf(nodeId)) {
            searchInFile(getLeafFile(nodeId), query, ansCount, result, resultLock, isTerminate);
            return;
        }
        searchInSubtree(leftChild(nodeId), query, ansCount, result, resultLock, isTerminate);
        searchInSubtree(rightChild(nodeId), query, ansCount, result, resultLock, isTerminate);
    }

    const std::filesystem::path &BallTree::getLeafFile(size_t nodeId) const {
        assert((1ull << _depth) - 1 <= nodeId);
        return _leafDataPaths[nodeId - (1ull << _depth) + 1];
    }

} // qtr