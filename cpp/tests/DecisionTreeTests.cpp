#include "DecisionTree.h"

#include "Common.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <queue>
#include <vector>
#include <utility>

using namespace qtr;

namespace {

class Less {
public:
    Less() : _n(0) {}
    Less(int n) : _n(n) {}
    bool operator()(int i) const { return i < _n; }
private:
    int _n;
};

} // anonymous namespace

TEST(DecisionTree, MAIN) {
    
    std::vector<int> data(10);
    for(int i = 0; i < data.size(); i++)
        data.at(i) = i;

    DecisionTree<::Less, std::vector<int>> tree;
    using Node = DecisionNode<::Less, std::vector<int>>;

    std::queue<Node *> nodes;

    tree.getRoot().setInfo(data);
    nodes.push(&tree.getRoot());

    while (!nodes.empty()) {

        Node *node = nodes.front();
        std::vector<int> &info = node->getInfo();

        if (info.size() >= 3) {

            std::size_t medianIdx = info.size() / 2;
            Node::Children children = node->setPred(::Less(info.at(medianIdx)));

            std::vector<int> infoTrue(info.begin(), info.begin() + medianIdx);
            children._true->setInfo(std::move(infoTrue));

            std::vector<int> infoFalse(info.begin() + medianIdx, info.end());
            children._false->setInfo(std::move(infoFalse));

            nodes.push({children._true});
            nodes.push({children._false});

            info.clear();
        }
        
        nodes.pop();
    }

    for(int i = 0; i < data.size(); i++) {
        const std::vector<int> &info = tree.search(i);
        auto it = std::find(info.begin(), info.end(), i);
        EXPECT_NE(it, info.end());
    }
}
