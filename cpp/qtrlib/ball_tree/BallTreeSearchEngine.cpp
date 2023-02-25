#include "BallTreeSearchEngine.h"

#include <cassert>
#include <numeric>
#include <random>
#include <future>
#include <chrono>

#include "fingerprint_table_io/FingerprintTableReader.h"
#include "Utils.h"

using namespace std;

namespace qtr {

    namespace {
        mt19937 random_generator(0);
    }

    void BallTreeSearchEngine::initLeafDataPaths() {
        size_t expectedFilesNumber = (_nodes.size() + 1) / 2;
        _leafDirPaths.resize(expectedFilesNumber);
        vector<bool> isInit(expectedFilesNumber, false);
        for (auto &dataDir: _dataDirectories) {
            for (auto &dirPath: filesystem::directory_iterator(dataDir)) {
                size_t nodeId = atoll(dirPath.path().filename().c_str());
                size_t index = nodeId - (1ull << _depth) + 1;
                assert(index < expectedFilesNumber);
                assert(!isInit[index]);
                isInit[index] = true;
                _leafDirPaths[index] = dirPath;
            }
        }
        assert(count(isInit.begin(), isInit.end(), false) == 0);
    }

    const filesystem::path &BallTreeSearchEngine::getLeafDir(size_t nodeId) const {
        assert((1ull << _depth) - 1 <= nodeId);
        return _leafDirPaths[nodeId - (1ull << _depth) + 1];
    }

    void BallTreeSearchEngine::processLeafGroup(BallTreeQueryData &queryData,
                                                vector <uint64_t> leafs, size_t group,
                                                size_t totalGroups) const {
        auto startTime = std::chrono::high_resolution_clock::now();
        queryData.tagStartTask();
        auto filterObject = queryData.getFilterObject();
        for (size_t i = group; i < leafs.size() && !queryData.checkShouldStop(); i += totalGroups) {
            auto res = searchInLeaf(leafs[i], queryData.getQueryFingerprint());
            if (res.empty())
                continue;
            filterObject->initBallTreeLeaf(getLeafDir(leafs[i]));
            queryData.filterAndAddAnswers(res, *filterObject);
        }
        std::chrono::duration<double> duration = std::chrono::high_resolution_clock::now() - startTime;
        ballTreeSearchTimer += duration.count();
        queryData.tagFinishTask();
    }

    void BallTreeSearchEngine::search(BallTreeQueryData &queryData, size_t threads) const {
        vector<uint64_t> leafs;
        auto startTime = std::chrono::high_resolution_clock::now();
        findLeafs(queryData.getQueryFingerprint(), root(), leafs);
        LOG(INFO) << "Search in " << leafs.size() << " leafs";
        shuffle(leafs.begin(), leafs.end(), random_generator);
        for (size_t i = 0; i < threads; i++) {
            auto task = std::async(std::launch::async, &BallTreeSearchEngine::processLeafGroup, this,
                                   std::ref(queryData), leafs, i, threads);
            queryData.addTask(std::move(task));
        }
        std::chrono::duration<double> duration = std::chrono::high_resolution_clock::now() - startTime;
        // multiply duration and threads to make time consumption percentage more accurate
        ballTreeSearchTimer += duration.count() * (double) threads;
    }

    vector <size_t> BallTreeSearchEngine::getLeafIds() const {
        vector<size_t> result((1ull << _depth));
        iota(result.begin(), result.end(), (1ull << _depth) - 1);
        return result;
    }

    void BallTreeSearchEngine::findLeafs(const IndigoFingerprint &fingerprint, size_t currentNode,
                                         vector <CIDType> &leafs) const {
        if (!(fingerprint <= _nodes[currentNode].centroid))
            return;
        if (isLeaf(currentNode)) {
            leafs.emplace_back(currentNode);
            return;
        }
        findLeafs(fingerprint, leftChild(currentNode), leafs);
        findLeafs(fingerprint, rightChild(currentNode), leafs);
    }

} // qtr