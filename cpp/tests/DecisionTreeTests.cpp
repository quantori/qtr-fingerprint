#include "DecisionTree.h"

#include "Common.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <queue>
#include <vector>
#include <utility>

using namespace qtr;

namespace {

class BitSet {
public:
    BitSet() : _pos(std::size_t(-1)) {}
    BitSet(std::size_t pos) : _pos(pos) {}
    bool operator()(uint8_t i) const { return (i >> _pos) & uint8_t(1); }
private:
    std::size_t _pos;
};

bool find(const std::vector<const std::vector<uint8_t> *> &infos, uint8_t query) {
    for(const std::vector<uint8_t> *info : infos) {
        if (std::find(info->begin(), info->end(), query) != info->end())
            return true;
    }
    return false;
}

} // anonymous namespace

TEST(DecisionTree, MAIN) {
    
    std::vector<uint8_t> data = {
        0b00001111, 0b00110011, 0b01010101, 0b10101010,
        0b11110000, 0b11001100, 0b10101010, 0b11100011,
        0b00011100, 0b11010101
    };

    DecisionTree<::BitSet, std::vector<uint8_t>> tree;
    using Node = DecisionNode<::BitSet, std::vector<uint8_t>>;

    std::queue<std::pair<std::size_t, Node *>> nodes;

    tree.getRoot().setInfo(data);
    nodes.push({0, &tree.getRoot()});

    while (!nodes.empty()) {

        size_t bit = nodes.front().first;
        Node *node = nodes.front().second;
        std::vector<uint8_t> &info = node->getInfo();

        if (bit < CHAR_BIT && info.size() >= 3) {

            ::BitSet pred(bit);
            
            Node::Children children = node->setPred(pred);
            auto it = std::partition(info.begin(), info.end(), pred);

            children._next->setInfo(std::vector<uint8_t>(info.begin(), it));
            children._false->setInfo(std::vector<uint8_t>(it, info.end()));

            nodes.push({bit + 1, children._next});
            nodes.push({bit + 1, children._false});

            info.clear();
        }
        
        nodes.pop();
    }

    for(uint8_t query : data) {
        auto infos = tree.search(query);
        EXPECT_TRUE(::find(infos, query));
    }
}
