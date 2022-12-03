#include "BallTreeSearchEngine.h"

#include <cassert>
#include <numeric>

#include "fingerprint_table_io/FingerprintTableReader.h"
#include "Utils.h"

using namespace std;

namespace qtr {

    void BallTreeSearchEngine::initLeafDataPaths() {
        size_t expectedFilesNumber = (_nodes.size() + 1) / 2;
        _leafDataPaths.resize(expectedFilesNumber);
        vector<bool> isInit(expectedFilesNumber, false);
        for (auto &dataDir: _dataDirectories) {
            for (auto &filePath: findFiles(dataDir, "")) {
                size_t nodeId = atoll(filePath.filename().c_str());
                size_t index = nodeId - (1ull << _depth) + 1;
                assert(index < expectedFilesNumber);
                assert(!isInit[index]);
                isInit[index] = true;
                _leafDataPaths[index] = filePath / ("data" + qtr::fingerprintTableExtension);
                assert(filesystem::is_regular_file(_leafDataPaths[index]));
            }
        }
        assert(count(isInit.begin(), isInit.end(), false) == 0);
    }

    const filesystem::path &BallTreeSearchEngine::getLeafFile(size_t nodeId) const {
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

    vector<future<void>>
    BallTreeSearchEngine::search(QueryData &queryData, size_t startDepth) const {
        vector<future<void>> tasks;
        for (size_t i = (1ull << (startDepth)) - 1; i < (1ull << (startDepth + 1)) - 1; i++) {
            tasks.emplace_back(
                    async(launch::async, &BallTreeSearchEngine::searchInSubtree, this, i,
                               ref(queryData))
            );
        }
        return tasks;
    }

    void BallTreeSearchEngine::putAnswer(CIDType ansValue, qtr::BallTreeSearchEngine::QueryData &queryData) {
        queryData.updateIsTerminate();
        if (queryData.isTerminate || !queryData.filter(ansValue)) {
            return;
        }
        queryData.addAnswer(ansValue);
    }

    vector<size_t> BallTreeSearchEngine::getLeafIds() const {
        vector<size_t> result((1ull << _depth));
        iota(result.begin(), result.end(), (1ull << _depth) - 1);
        return result;
    }

    void BallTreeSearchEngine::QueryData::updateIsTerminate() {
        isTerminate |= result.size() >= ansCount;
    }

    void BallTreeSearchEngine::QueryData::addAnswer(CIDType value) {
        if (result.size() < ansCount) {
            lock_guard<mutex> lock(resultLock);
            if (result.size() < ansCount) {
                result.emplace_back(value);
            }
        }
    }

    BallTreeSearchEngine::QueryData::QueryData(size_t ansCount,
                                               IndigoFingerprint query,
                                               function<bool(CIDType)> filter) :
            ansCount(ansCount),
            query(query),
            filter(move(filter)),
            isTerminate(false) {
    }
} // qtr