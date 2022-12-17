#include "BallTreeSearchEngine.h"

#include <cassert>
#include <numeric>
#include <random>
#include <future>

#include "fingerprint_table_io/FingerprintTableReader.h"
#include "Utils.h"

using namespace std;

namespace qtr {

    namespace {
        mt19937 random_generator(0);
    }

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

    void BallTreeSearchEngine::processLeafGroup(BallTreeQueryData &queryData,
                                                vector <uint64_t> leafs, size_t group,
                                                size_t totalGroups) const {
        for (size_t i = group; i < leafs.size(); i += totalGroups) {
            auto res = searchInLeaf(leafs[i], queryData.getQueryFingerprint());
            queryData.filterAndAddAnswers(res);
        }
    }

    void BallTreeSearchEngine::search(BallTreeQueryData &queryData, size_t threads) const {
        vector<uint64_t> leafs;
        findLeafs(queryData.getQueryFingerprint(), root(), leafs);
        LOG(INFO) << "Search in " << leafs.size() << " leafs";
        shuffle(leafs.begin(), leafs.end(), random_generator);
        for (size_t i = 0; i < threads; i++) {
            auto task = std::async(std::launch::async, &BallTreeSearchEngine::processLeafGroup, this,
                                   std::ref(queryData), leafs, i, threads);
            queryData.addTask(std::move(task));
        }
    }

    vector <size_t> BallTreeSearchEngine::getLeafIds() const {
        vector<size_t> result((1ull << _depth));
        iota(result.begin(), result.end(), (1ull << _depth) - 1);
        return result;
    }

    void
    BallTreeSearchEngine::findLeafs(const IndigoFingerprint &fingerprint, size_t currentNode,
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