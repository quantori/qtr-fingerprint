#include "BallTreeSearchEngine.h"

#include <cassert>
#include <numeric>
#include <random>
#include <future>
#include <chrono>

#include "fingerprint_table_io/FingerprintTableReader.h"
#include "Utils.h"
#include "QtrBallTreeLeafInitMixin.h"
#include "Profiling.h"

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
    BallTreeSearchEngine::processLeafGroup(QueryDataWithFingerprint &queryData, const vector <uint64_t> &leaves) const {
        ProfileScope("processLeafGroup");
        auto filterObject = queryData.getFilterObject();
        for (size_t i = 0; i < leaves.size() && !queryData.checkShouldStop(); i++) {
            processLeaf(queryData, leaves[i], *filterObject);
        }
        if (queryData.checkTimeOut()) {
            LOG(INFO) << "Search stopped due to timeout";
        }
        queryData.tagFinishTask();
    }

    void BallTreeSearchEngine::search(QueryDataWithFingerprint &queryData, size_t threads) const {
        if (threads == 1) {
            queryData.addTask(async(launch::async, [this, &queryData] {
                auto filterObject = queryData.getFilterObject();
                oneThreadSearch(queryData, root(), *filterObject);
            }));
        } else {
            ProfileScope("Search fingerprints in BallTree");
            vector<uint64_t> leaves;
            findLeaves(queryData.getQueryFingerprint(), root(), leaves);
            LOG(INFO) << "Search in " << leaves.size() << " leaves";
            auto leafGroups = divideLeavesIntoGroups(leaves, threads);
            for (auto &leafGroup: leafGroups) {
                auto task = async(launch::async, &BallTreeSearchEngine::processLeafGroup, this,
                                  ref(queryData), leafGroup);
                queryData.addTask(std::move(task));
            }
        }
    }

    vector <size_t> BallTreeSearchEngine::getLeafIds() const {
        vector<size_t> result((1ull << _depth));
        iota(result.begin(), result.end(), (1ull << _depth) - 1);
        return result;
    }

    void BallTreeSearchEngine::findLeaves(const Fingerprint &fingerprint, size_t currentNode,
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

    void BallTreeSearchEngine::oneThreadSearch(QueryDataWithFingerprint &queryData, size_t currentNode,
                                               ByIdAnswerFilter &filterObject) const {
        if (queryData.checkFoundEnoughAnswers()) {
            return;
        }
        if (!(queryData.getQueryFingerprint() <= _nodes[currentNode].centroid))
            return;
        if (isLeaf(currentNode)) {
            if ((currentNode & ((1 << 10) - 1)) == 0 &&
                (queryData.checkTimeOut() || queryData.checkShouldStop())) { // check timeout once per 1024 leaves
                return;
            }
            processLeaf(queryData, currentNode, filterObject);
            return;
        }
        oneThreadSearch(queryData, leftChild(currentNode), filterObject);
        oneThreadSearch(queryData, rightChild(currentNode), filterObject);
    }

    void BallTreeSearchEngine::processLeaf(QueryDataWithFingerprint &queryData, uint64_t leafId,
                                           ByIdAnswerFilter &filterObject) const {
        auto res = searchInLeaf(leafId, queryData.getQueryFingerprint());
        if (res.empty())
            return;
        tryInitBallTreeLeaf(filterObject, getLeafDir(leafId));
        queryData.filterAndAddAnswers(res, filterObject);
    }

} // qtr