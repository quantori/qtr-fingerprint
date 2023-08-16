#include "BallTreeSearchEngine.h"

#include <cassert>
#include <numeric>
#include <random>
#include <future>
#include <chrono>
#include <algorithm>

#include "fingerprint_table_io/FingerprintTableReader.h"
#include "Utils.h"
#include "QtrBallTreeLeafInitMixin.h"

using namespace std;

namespace qtr {

    void BallTreeSearchEngine::initLeafDataPaths() {
        size_t expectedFilesNumber = (_nodes.size() + 1) / 2;
        _leafDirPaths.resize(expectedFilesNumber);
        vector<bool> isInit(expectedFilesNumber, false);
        for (auto &dataDir: _dataDirectories) {
            for (auto &dirPath: filesystem::directory_iterator(dataDir)) {
                size_t nodeId = stoll(dirPath.path().filename().string());
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

    void
    BallTreeSearchEngine::processLeafGroup(BallTreeQueryData &queryData, const vector <uint64_t> &leaves) const {
        auto startTime = chrono::high_resolution_clock::now();
        auto filterObject = queryData.getFilterObject();
        for (size_t i = 0; i < leaves.size() && !queryData.checkShouldStop(); i++) {
            auto res = searchInLeaf(leaves[i], queryData.getQueryFingerprint());
            if (res.empty())
                continue;
            tryInitBallTreeLeaf(*filterObject, getLeafDir(leaves[i]));
            queryData.filterAndAddAnswers(res, *filterObject);
        }
        if (queryData.checkTimeOut()) {
            LOG(INFO) << "Search stopped due to timeout";
        }

        chrono::duration<double> duration = chrono::high_resolution_clock::now() - startTime;
        ballTreeSearchTimer += duration.count();
        queryData.tagFinishTask();
    }

    void BallTreeSearchEngine::search(BallTreeQueryData &queryData, size_t threads) const {
        vector<uint64_t> leaves;
        auto startTime = chrono::high_resolution_clock::now();
        findLeaves(queryData.getQueryFingerprint(), root(), leaves);
        LOG(INFO) << "Search in " << leaves.size() << " leaves";
        auto leafGroups = divideLeavesIntoGroups(leaves, threads);
        for (auto &leafGroup: leafGroups) {
            auto task = async(launch::async, &BallTreeSearchEngine::processLeafGroup, this,
                              ref(queryData), leafGroup);
            queryData.addTask(std::move(task));
        }
        chrono::duration<double> duration = chrono::high_resolution_clock::now() - startTime;
        // multiply duration and groups number to make time consumption percentage more accurate
        ballTreeSearchTimer += duration.count() * (double) leafGroups.size();
    }

    vector <size_t> BallTreeSearchEngine::getLeafIds() const {
        vector<size_t> result((1ull << _depth));
        iota(result.begin(), result.end(), (1ull << _depth) - 1);
        return result;
    }

    void BallTreeSearchEngine::findLeaves(const IndigoFingerprint &fingerprint, size_t currentNode,
                                          vector <CIDType> &leaves) const {
        if (!(fingerprint <= _nodes[currentNode].centroid))
            return;
        if (isLeaf(currentNode)) {
            leaves.emplace_back(currentNode);
            return;
        }
        findLeaves(fingerprint, leftChild(currentNode), leaves);
        findLeaves(fingerprint, rightChild(currentNode), leaves);
    }

} // qtr