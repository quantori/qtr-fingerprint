#include <utility>
#include <future>
#include <cstdlib>
#include <random>

#include "BallTreeBuilder.h"
#include "data_io/fingerprint_table_io/FingerprintTableReader.h"
#include "data_io/fingerprint_table_io/FingerprintTableWriter.h"

namespace qtr {

    namespace {

        const std::string dataFilename = "data.ft";
        std::mt19937 randomGenerator(0);

        void splitFileByBit(const std::filesystem::path &file, size_t splitBit, const std::filesystem::path &leftFile,
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
            LOG(INFO) << "Delete " << file;
            std::filesystem::remove(file);
        }

        void prepareManyDirsForNode(const std::vector<std::filesystem::path> &dataDirectories, size_t nodeId) {
            for (auto &dirPath: dataDirectories) {
                LOG(INFO) << "Create directory " << dirPath / std::to_string(nodeId);
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
                LOG(INFO) << "Delete " << source;
                std::filesystem::remove(source);
            }
        }

    } // namespace

    BallTreeBuilder::BallTreeBuilder(size_t depth, size_t parallelizationDepth,
                                     std::vector<std::filesystem::path> dataDirectories,
                                     const BitSelector &bitSelector) : BallTree(std::move(dataDirectories)) {
        _depth = depth;
        _nodes.resize((1ull << (depth + 1)) - 1);
        buildTree(depth, parallelizationDepth, bitSelector);
    }

    void BallTreeBuilder::buildTree(size_t depth, size_t parallelizationDepth, const BitSelector &bitSelector) {
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
            for (auto &dir: _dataDirectories) {
                auto nodeDir = dir / std::to_string(nodes[i]);
                if (nodeDir != destinationPath.parent_path() &&
                    std::filesystem::exists(nodeDir)) {
                    LOG(INFO) << "Delete " << nodeDir;
                    std::filesystem::remove_all(nodeDir);
                }
            }
        }
        LOG(INFO) << "Finish subtree parallelization preparation";
        LOG(INFO) << "Start last levels building";
        std::vector<std::future<void>> tasks;
        for (size_t i = 0; i < nodes.size(); i++) {
            tasks.emplace_back(std::async(std::launch::async, &BallTreeBuilder::buildLastLevels, this, nodes[i],
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

    void BallTreeBuilder::splitNodeManyFiles(size_t nodeId, const BitSelector &bitSelector) {
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
                    std::async(std::launch::async, splitFileByBit, nodeFiles[i], splitBit, leftFile, rightFile,
                               splitPolicy[i].first, splitPolicy[i].second)
            );
        }
        for (auto &task: tasks) {
            task.get();
        }
        LOG(INFO) << "Finish parallel splitting of node " << nodeId;
        deleteNodeFromFilesystem(nodeId);
    }

    void BallTreeBuilder::splitNodeOneFile(size_t nodeId, const BitSelector &bitSelector,
                                           const std::filesystem::path &dataDirectory) {
        std::filesystem::path nodeFile = dataDirectory / std::to_string(nodeId) / dataFilename;
        auto [_, stat] = collectFilesStatistic({nodeFile});
        size_t splitBit = bitSelector(stat);
        auto [leftSize, rightSize] = findSplitPolicy(stat, splitBit);
        std::filesystem::path leftFile = dataDirectory / std::to_string(leftChild(nodeId)) / dataFilename;
        std::filesystem::path rightFile = dataDirectory / std::to_string(rightChild(nodeId)) / dataFilename;
        std::filesystem::create_directory(leftFile.parent_path());
        std::filesystem::create_directory(rightFile.parent_path());
        splitFileByBit(nodeFile, splitBit, leftFile, rightFile, leftSize, rightSize);
        LOG(INFO) << "Delete " << nodeFile.parent_path();
        std::filesystem::remove_all(nodeFile.parent_path());
    }


    std::vector<size_t>
    BallTreeBuilder::buildFirstLevels(size_t nodeId, size_t levels, const BitSelector &bitSelector) {
        if (levels == 0)
            return {nodeId};
        splitNodeManyFiles(nodeId, bitSelector);
        std::vector<size_t> nodes = buildFirstLevels(leftChild(nodeId), levels - 1, bitSelector);
        std::vector<size_t> rightNodes = buildFirstLevels(rightChild(nodeId), levels - 1, bitSelector);
        nodes.insert(nodes.end(), rightNodes.begin(), rightNodes.end());
        return nodes;
    }

    std::vector<std::filesystem::path> BallTreeBuilder::getNodeFiles(size_t nodeId) const {
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

    void BallTreeBuilder::buildLastLevels(size_t nodeId, size_t levels, const BitSelector &bitSelector,
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

    void BallTreeBuilder::calculateCentroid(size_t nodeId) {
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

    void BallTreeBuilder::deleteNodeFromFilesystem(size_t nodeId) const {
        for (auto &dir: _dataDirectories) {
            auto nodeDir = dir / std::to_string(nodeId);
            if (std::filesystem::exists(nodeDir)) {
                LOG(INFO) << "Delete " << nodeDir;
                std::filesystem::remove_all(nodeDir);
            }
        }
    }

} // qtr