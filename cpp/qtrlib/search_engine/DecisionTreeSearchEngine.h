#pragma once

#include "DecisionTree.h"
#include "FingerprintTable.h"
#include "SearchEngineInterface.h"

#include "IndigoMolecule.h"

#include <vector>

namespace qtr {

class DecisionTreeSearchEngine : public SearchEngineInterface {
public:
    DecisionTreeSearchEngine() = delete;
    explicit DecisionTreeSearchEngine(const indigo_cpp::IndigoSessionPtr &indigoSessionPtr);

    ~DecisionTreeSearchEngine() override;

    void build(const std::string &path) override;

    std::vector<indigo_cpp::IndigoMolecule> findOverMolecules(const indigo_cpp::IndigoQueryMolecule &mol) override;

private:
    static const size_t maxLeafSize = 100;

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

    DecisionTree<BitSet, IndigoFingerprintTableView> _decisionTree;
    IndigoFingerprintTable _fingerprintTable;
    std::vector<indigo_cpp::IndigoMolecule> _molecules;
    indigo_cpp::IndigoSessionPtr _indigoSessionPtr;
};

} // namespace qtr