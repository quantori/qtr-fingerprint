#pragma once

#include <vector>
#include <memory>
#include <cmath>
#include <cassert>
#include <numeric>

#include "frameworks/FrameworkInterface.h"
#include "CachedDataset.h"

template<typename FrameworkT> requires FrameworkInterface<FrameworkT>
class BallTreeSplitter {
public:
    BallTreeSplitter(CachedDataset<FrameworkT> &dataset, size_t rootId, size_t nodesCount) :
            _nodeIndices(nodesCount),
            _dataset(dataset),
            _splitBits(nodesCount) {
        auto &initIndices = _nodeIndices.at(rootId);
        initIndices.resize(_dataset.size());
        std::iota(initIndices.begin(), initIndices.end(), size_t(0));
    }

    void split(size_t sourceNode, size_t leftNode, size_t rightNode) {
        size_t splitBit = findSplitBit(sourceNode);
        _splitBits[sourceNode] = splitBit;
        splitByBit(sourceNode, leftNode, rightNode, splitBit);
    }

    std::vector<size_t> &&extractIndices(size_t nodeId) {
        return std::move(_nodeIndices.at(nodeId));
    }

    size_t getSplitBit(size_t nodeId) {
        return _splitBits.at(nodeId);
    }

private:
    std::vector<std::vector<size_t>> _nodeIndices;
    CachedDataset<FrameworkT> &_dataset;
    std::vector<size_t> _splitBits;

    size_t findSplitBit(size_t nodeId) {
        auto hist = getPerBitHistogram(nodeId);
        size_t splitBit = selectSplitBitByHistogram(hist);
        return splitBit;
    }

    void splitByBit(size_t sourceNode, size_t leftNode, size_t rightNode, size_t splitBit) {
        auto &sourceIndices = _nodeIndices.at(sourceNode);
        auto &leftIndices = _nodeIndices.at(leftNode);
        auto &rightIndices = _nodeIndices.at(rightNode);
        size_t n = sourceIndices.size();
        for (auto &idx: sourceIndices) {
            if (leftIndices.size() * 2 >= n) {
                rightIndices.push_back(idx);
                continue;
            } else if (rightIndices.size() * 2 >= n) {
                leftIndices.push_back(idx);
                continue;
            }
            auto &fingerprint = _dataset.fingerprint(idx);
            bool splitBitVal = FrameworkT::getFingerprintBit(fingerprint, splitBit);
            if (splitBitVal) {
                rightIndices.push_back(idx);
            } else {
                leftIndices.push_back(idx);
            }
        }
        auto lSize = leftIndices.size();
        auto rSize = rightIndices.size();
        assert(lSize + rSize == n);
        assert(abs((int) lSize - (int) rSize) <= 1);
        sourceIndices.clear();
        sourceIndices.shrink_to_fit();
//        LOG(INFO) << "Split node " << sourceNode << "(" << n << ") into " << leftNode << "(" << lSize << ") and "
//                  << rightNode << "(" << rSize << ")";
    }

    std::vector<size_t> getPerBitHistogram(size_t nodeId) {
        std::vector<size_t> hist(FrameworkT::getFingerprintSize());
        for (size_t idx: _nodeIndices.at(nodeId)) {
            const auto &fingerprint = _dataset.fingerprint(idx);
            for (size_t j = 0; j < FrameworkT::getFingerprintSize(); j++) {
                hist[j] += FrameworkT::getFingerprintBit(fingerprint, j);
            }
        }
        return hist;
    }

    size_t selectSplitBitByHistogram(const std::vector<size_t> &hist) {
        size_t splitBit = 0;
        size_t bestBalance = getBalanceAtBit(hist, splitBit);
        for (size_t i = 1; i < hist.size(); i++) {
            auto currentBalance = getBalanceAtBit(hist, i);
            if (bestBalance > currentBalance) {
                bestBalance = currentBalance;
                splitBit = i;
            }
        }
        return splitBit;
    }

    size_t getBalanceAtBit(const std::vector<size_t> &hist, size_t idx) {
        auto ones = (int) hist.at(idx);
        auto zeros = (int) _dataset.size() - ones;
        return std::abs(ones - zeros);
    }
};
