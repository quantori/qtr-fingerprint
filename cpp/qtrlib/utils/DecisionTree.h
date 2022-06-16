#pragma once

#include "FingerprintTableView.h"

#include <memory>
#include <utility>

namespace qtr {

template<class Pred, class Info>
class DecisionNode {
public:
    using Node = DecisionNode<Pred, Info>;

    struct Children {
        Node *_true;
        Node *_false;
    };

    Children setPred(const Pred &pred) {
        _pred = pred;
        _true = std::make_unique<Node>();
        _false = std::make_unique<Node>();
        return {_true.get(), _false.get()};
    }

    void setInfo(Info &&info) {
        _info = std::move(info);
    }

    void setInfo(const Info &info) {
        _info = info;
    }

    bool isLeaf() const {
        return !_true;
    }

    Info &getInfo() {
        return _info;
    }

    const Info &getInfo() const {
        return _info;
    }

    template<class Query>
    const Node *next(const Query &query) const {
        if (_pred(query))
            return _true.get();
        return _false.get();
    }

private:
    Pred _pred;
    Info _info;

    std::unique_ptr<Node> _true;
    std::unique_ptr<Node> _false;
};

template<class Pred, class Info>
class DecisionTree {
public:
    using Node = DecisionNode<Pred, Info>;

    Node &getRoot() { return _root; }

    template<class Query>
    const Info &search(const Query &query) const {
        const Node *node = &_root;
        while(!node->isLeaf())
            node = node->next(query);
        return node->getInfo();
    }

private:  
    Node _root;
};

} // namespace qtr