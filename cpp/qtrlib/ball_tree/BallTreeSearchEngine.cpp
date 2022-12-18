#include "BallTreeSearchEngine.h"

#include <cassert>
#include <numeric>
#include <random>

#include "fingerprint_table_io/FingerprintTableReader.h"
#include "Utils.h"

namespace qtr {

    namespace {
        std::mt19937 random_generator(0);
    }

    void BallTreeSearchEngine::initLeafDataPaths() {
        size_t expectedFilesNumber = (_nodes.size() + 1) / 2;
        _leafDataPaths.resize(expectedFilesNumber);
        std::vector<bool> isInit(expectedFilesNumber, false);
        for (auto &dataDir: _dataDirectories) {
            for (auto &filePath: findFiles(dataDir, "")) {
                size_t nodeId = atoll(filePath.filename().c_str());
                size_t index = nodeId - (1ull << _depth) + 1;
                assert(index < expectedFilesNumber);
                assert(!isInit[index]);
                isInit[index] = true;
                _leafDataPaths[index] = filePath / ("data" + qtr::fingerprintTableExtension);
                assert(std::filesystem::is_regular_file(_leafDataPaths[index]));
            }
        }
        assert(std::count(isInit.begin(), isInit.end(), false) == 0);
    }

    const std::filesystem::path &BallTreeSearchEngine::getLeafFile(size_t nodeId) const {
        assert((1ull << _depth) - 1 <= nodeId);
        return _leafDataPaths[nodeId - (1ull << _depth) + 1];
    }

    void BallTreeSearchEngine::searchInSubtree(size_t nodeId, QueryData &queryData) const {
        if (queryData.isTerminate || !(queryData.query <= _nodes[nodeId].centroid)) {
            return;
        }
        if (isLeaf(nodeId)) {
            searchInLeaf(nodeId, queryData);
            return;
        }
        searchInSubtree(leftChild(nodeId), queryData);
        searchInSubtree(rightChild(nodeId), queryData);
    }

    void
    BallTreeSearchEngine::findLeafs(const IndigoFingerprint &fingerprint, size_t currentNode,
                                    std::vector <CIDType> &leafs) const {
        if (!(fingerprint <= _nodes[currentNode].centroid))
            return;
        if (isLeaf(currentNode)) {
            leafs.emplace_back(currentNode);
            return;
        }
        findLeafs(fingerprint, leftChild(currentNode), leafs);
        findLeafs(fingerprint, rightChild(currentNode), leafs);
    }

    void BallTreeSearchEngine::processLeafGroup(QueryData &queryData,
                                                std::vector <uint64_t> leafs, size_t group,
                                                size_t totalGroups) const {
        for (size_t i = group; i < leafs.size(); i += totalGroups) {
            if (queryData.isTerminate)
                break;
            searchInLeaf(leafs[i], queryData);
        }
    }

    std::vector<CIDType>
    BallTreeSearchEngine::search(const IndigoFingerprint &query, size_t ansCount, size_t startDepth,
                                 const std::function<bool(CIDType)> &filter) const {

        std::vector<uint64_t> leafs;
        findLeafs(query, root(), leafs);
        LOG(INFO) << "Search in " << leafs.size() << " leafs";
        shuffle(leafs.begin(), leafs.end(), random_generator);

        std::vector<std::future<void>> tasks;
        std::vector<CIDType> result;
        std::mutex resultsLock;
        const size_t threads = 16;
        QueryData queryData = {query, result, resultsLock, ansCount, false, filter};
        for (size_t i = 0; i < threads; i++) {
            tasks.emplace_back(std::async(std::launch::async, &BallTreeSearchEngine::processLeafGroup, this,
                                          std::ref(queryData), leafs, i, threads));
        }
        for (auto &task: tasks) {
            task.get();
        }
        return result;
    }

    void BallTreeSearchEngine::putAnswer(CIDType ansValue, qtr::BallTreeSearchEngine::QueryData &queryData) {
        queryData.updateIsTerminate();
        if (queryData.isTerminate || !queryData.filter(ansValue)) {
            return;
        }
        queryData.addAnswer(ansValue);
    }

    std::vector<size_t> BallTreeSearchEngine::getLeafIds() const {
        std::vector<size_t> result((1ull << _depth));
        std::iota(result.begin(), result.end(), (1ull << _depth) - 1);
        return result;
    }

    void BallTreeSearchEngine::QueryData::updateIsTerminate() {
        isTerminate |= result.size() >= ansCount;
    }

    void BallTreeSearchEngine::QueryData::addAnswer(CIDType value) {
        std::lock_guard<std::mutex> lock(resultLock);
        result.emplace_back(value);
    }
} // qtr