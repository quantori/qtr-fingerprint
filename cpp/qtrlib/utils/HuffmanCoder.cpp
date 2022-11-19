#include "HuffmanCoder.h"

#include <glog/logging.h>

#include <queue>
#include <functional>
#include <cassert>

namespace qtr {

    HuffmanCoder::HuffmanCoder(const std::map<char, uint64_t>& symbolsFrequency) {
        createTree(symbolsFrequency);
        createSymbolsMap();
    }

    void HuffmanCoder::createTree(const std::map<char, uint64_t>& symbolsFrequency) {
        std::priority_queue<std::pair<uint64_t, size_t>,
                std::vector<std::pair<uint64_t, size_t>>,
                std::greater<>> freeNodes;
        assert(!symbolsFrequency.empty() && "can't build tree for empty alphabet");
        for (auto& [symbol, priority] : symbolsFrequency) {
            freeNodes.emplace(priority, _treeNodes.size());
            _treeNodes.emplace_back(symbol);
        }
        assert(!freeNodes.empty());
        while (freeNodes.size() > 1) {
            auto [priority1, id1] = freeNodes.top();
            freeNodes.pop();
            auto [priority2, id2] = freeNodes.top();
            freeNodes.pop();
            freeNodes.emplace(priority1 + priority2, _treeNodes.size());
            _treeNodes.emplace_back('\0', id1, id2);
        }
        _treeRoot = freeNodes.top().second;
    }

    void HuffmanCoder::createSymbolsMap() {
        std::function<void(size_t, std::vector<bool> &)> DFS;
        DFS = [this, &DFS](size_t nodeId, std::vector<bool> &code) {
            if (_treeNodes[nodeId].symbol != '\0') {
                getSymbolCode(_treeNodes[nodeId].symbol) = code;
                getSymbolCode(_treeNodes[nodeId].symbol).shrink_to_fit();
                return;
            }
            code.push_back(false);
            DFS(_treeNodes[nodeId].left, code);
            code.back() = true;
            DFS(_treeNodes[nodeId].right, code);
            code.pop_back();
        };
        std::fill(std::begin(_symbolsMap), std::end(_symbolsMap), std::vector<bool>());
        std::vector<bool> code;
        DFS(_treeRoot, code);
    }

    std::vector<bool> HuffmanCoder::encode(const std::string &string) const {
        std::vector<bool> result;
        for (auto &symbol: string) {
            auto &symbolCode = getSymbolCode(symbol);
            assert(!symbolCode.empty());
            result.insert(result.end(), symbolCode.begin(), symbolCode.end());
        }
        result.shrink_to_fit();
        return result;
    }

    std::string HuffmanCoder::decode(const std::vector<bool> &code) const {
        std::string result;
        size_t nodeId = _treeRoot;
        for (bool b: code) {
            auto &currentNode = _treeNodes[nodeId];
            nodeId = b ? currentNode.right : currentNode.left;
            auto &nextNode = _treeNodes[nodeId];
            if (nextNode.symbol == '\0')
                continue;
            result += nextNode.symbol;
            nodeId = _treeRoot;
        }
        assert(nodeId == _treeRoot);
        return result;
    }

    const std::vector<bool> &HuffmanCoder::getSymbolCode(char symbol) const {
        uint8_t unsignedSymbol = symbol;
        return _symbolsMap[unsignedSymbol];
    }

    std::vector<bool> &HuffmanCoder::getSymbolCode(char symbol) {
        return const_cast<std::vector<bool> &>(
                const_cast<const HuffmanCoder *>(this)->getSymbolCode(symbol)
        );
    }

    HuffmanCoder::Builder &HuffmanCoder::Builder::operator+=(const std::string &s) {
        for (char symbol: s) {
            symbolsFrequency[symbol]++;
        }
        return *this;
    }

    HuffmanCoder::Builder &HuffmanCoder::Builder::operator+=(const HuffmanCoder::Builder &other) {
        for (auto &[symbol, frequency]: other.symbolsFrequency) {
            symbolsFrequency[symbol] += frequency;
        }
        return *this;
    }

    HuffmanCoder HuffmanCoder::Builder::build() const {
        return { symbolsFrequency };
    }
} // qtr