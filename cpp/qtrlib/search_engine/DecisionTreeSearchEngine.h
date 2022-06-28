#pragma once

#include "DecisionTree.h"
#include "FingerprintTable.h"
#include "FingerprintTableSearchEngine.h"

#include "IndigoMolecule.h"

#include <vector>

namespace qtr {

template<class SplittingStrategy>
class DecisionTreeSearchEngine : public FingerprintTableSearchEngine {
public:
    DecisionTreeSearchEngine(
        const indigo_cpp::IndigoSessionPtr &indigoSessionPtr,
        size_t maxLeafSize = 100);

    void build(const std::string &path) override;

private:
    std::vector<const IndigoFingerprintTableView *> findTableViews(const qtr::IndigoFingerprint &fp) const override;

    class BitSet {
    public:
        BitSet() : _pos(std::size_t(-1)) {}
        BitSet(std::size_t pos) : _pos(pos) {}

        bool operator ()(const IndigoFingerprint &fp) const {
            return fp.test(_pos);
        };
    private:
        std::size_t _pos;
    };

    const size_t _maxLeafSize;
    SplittingStrategy _splittingStrategy;
    DecisionTree<BitSet, IndigoFingerprintTableView> _decisionTree;
};

class SplittingStrategyTrivial {
public:
    std::size_t operator()(std::size_t bitIndex, const IndigoFingerprintTableView &view) {
        return bitIndex;
    }
};

class SplittingStrategyOptimal {
public:
    std::size_t operator()(std::size_t bitIndex, const IndigoFingerprintTableView &view);
};

extern template class DecisionTreeSearchEngine<SplittingStrategyTrivial>;
extern template class DecisionTreeSearchEngine<SplittingStrategyOptimal>;

} // namespace qtr