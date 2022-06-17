#pragma once

#include "DecisionTree.h"
#include "FingerprintTable.h"
#include "SearchEngineInterface.h"

#include "IndigoMolecule.h"

#include <vector>

namespace qtr {

template<class SplittingStrategy>
class DecisionTreeSearchEngine : public SearchEngineInterface {
public:
    DecisionTreeSearchEngine() = delete;
    
    explicit DecisionTreeSearchEngine(
        const indigo_cpp::IndigoSessionPtr &indigoSessionPtr,
        size_t maxLeafSize = 100);

    ~DecisionTreeSearchEngine() override;

    void build(const std::string &path) override;

    std::vector<indigo_cpp::IndigoMolecule> findOverMolecules(const indigo_cpp::IndigoQueryMolecule &mol) override;

private:
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
    DecisionTree<BitSet, IndigoFingerprintTableView> _decisionTree;
    IndigoFingerprintTable _fingerprintTable;
    std::vector<indigo_cpp::IndigoMolecule> _molecules;
    indigo_cpp::IndigoSessionPtr _indigoSessionPtr;
    SplittingStrategy _splittingStrategy;
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