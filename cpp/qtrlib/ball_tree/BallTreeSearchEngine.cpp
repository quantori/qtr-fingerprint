#include "BallTreeSearchEngine.h"

#include <cassert>

#include "fingerprint_table_io/FingerprintTableReader.h"
#include "Utils.h"

namespace qtr {

    namespace {
        void putAnswer(std::vector<size_t> &result, size_t value, size_t ansCount, std::mutex &resultLock,
                       bool &isTerminate, const std::function<bool(size_t)>& filter) {
            isTerminate |= result.size() >= ansCount;
            if (isTerminate || !filter(value)) {
                return;
            }
            std::lock_guard<std::mutex> lock(resultLock);
            if (result.size() < ansCount) {
                result.emplace_back(value);
            }
        }

        void searchInFile(const std::filesystem::path &filePath, const IndigoFingerprint &query, size_t ansCount,
                          std::vector<size_t> &result, std::mutex &resultsLock, bool &isTerminate,
                          const std::function<bool(size_t)> &filter) {
            FingerprintTableReader reader(filePath);
            for (const auto &[id, fingerprint]: reader) {
                if (query <= fingerprint)
                    putAnswer(result, id, ansCount, resultsLock, isTerminate, filter);
            }
        }
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

    void BallTreeSearchEngine::searchInSubtree(size_t nodeId, const IndigoFingerprint &query, size_t ansCount,
                                               std::vector<size_t> &result,
                                               std::mutex &resultLock, bool &isTerminate,
                                               const std::function<bool(size_t)> &filter) const {
        if (isTerminate || !(query <= _nodes[nodeId].centroid)) {
            return;
        }
        if (isLeaf(nodeId)) {
            searchInFile(getLeafFile(nodeId), query, ansCount, result, resultLock, isTerminate, filter);
            return;
        }
        searchInSubtree(leftChild(nodeId), query, ansCount, result, resultLock, isTerminate, filter);
        searchInSubtree(rightChild(nodeId), query, ansCount, result, resultLock, isTerminate, filter);
    }

    std::vector<size_t>
    BallTreeSearchEngine::search(const IndigoFingerprint &query, size_t ansCount, size_t startDepth,
                                 const std::function<bool(size_t)> &filter) const {
        std::vector<std::future<void>> tasks;
        bool isTerminate = false;
        std::vector<size_t> results;
        std::mutex resultsLock;
        for (size_t i = (1ull << (startDepth)) - 1; i < (1ull << (startDepth + 1)) - 1; i++) {
            tasks.emplace_back(
                    std::async(std::launch::async, &BallTreeSearchEngine::searchInSubtree, this, i, std::cref(query),
                               ansCount,
                               std::ref(results), std::ref(resultsLock), std::ref(isTerminate), std::ref(filter))
            );
        }
        for (auto &task: tasks) {
            task.get();
        }
        return results;
    }

} // qtr