#include "DecisionTreeSearchEngine.h"

#include "Histogram.h"
#include "QtrIndigoFingerprint.h"

#include "IndigoSDFileIterator.h"
#include "IndigoSubstructureMatcher.h"

#include "indigo.h"

#include <glog/logging.h>

#include <algorithm>
#include <queue>
#include <utility>

using namespace indigo_cpp;

namespace qtr{

template<class SplittingStrategy>
DecisionTreeSearchEngine<SplittingStrategy>::DecisionTreeSearchEngine(
    const IndigoSessionPtr &indigoSessionPtr, size_t maxLeafSize)
    : _indigoSessionPtr(indigoSessionPtr)
    , _maxLeafSize(maxLeafSize)
{}

template<class SplittingStrategy>
DecisionTreeSearchEngine<SplittingStrategy>::~DecisionTreeSearchEngine()
{}

template<class SplittingStrategy>
void DecisionTreeSearchEngine<SplittingStrategy>::build(const std::string &path)
{
    _molecules.clear();
    _fingerprintTable.clear();

    size_t moleculesNumber = 0;
    IndigoSDFileIterator iterator = _indigoSessionPtr->iterateSDFile(path);

    for(IndigoMoleculeSPtr &molecule : iterator) {
        
        molecule->aromatize();
        _molecules.push_back(std::move(*molecule));

        QtrIndigoFingerprint fingerprint(_molecules.back(), "sub");
        _fingerprintTable.push_back(qtr::IndigoFingerprint());
        _fingerprintTable.back().setBytes(fingerprint.data());

        moleculesNumber++;
        if (moleculesNumber % 1000 == 0)
            LOG(INFO) << "Processed " << moleculesNumber << " molecules...";
    }

    ////////////////////////

    // Histogram histogram(CHAR_BIT*qtr::IndigoFingerprint::sizeInBytes);
    // for(const qtr::IndigoFingerprint &fp : _fingerprintTable) {
    //     for(size_t bit = 0; bit < fp.size(); bit++)
    //         histogram.add(bit, Histogram::CounterType(fp.test(bit)));
    // }

    // std::vector<size_t> bitsPerm(CHAR_BIT*qtr::IndigoFingerprint::sizeInBytes);
    // for(size_t i = 0; i < bitsPerm.size(); i++)
    //     bitsPerm.at(i) = i;

    // std::sort(bitsPerm.begin(), bitsPerm.end(), [&histogram](size_t left, size_t right) {
    //     return histogram.bins().at(left) < histogram.bins().at(right);
    // });

    // for(size_t i = 0; i < bitsPerm.size(); i++) {
    //     size_t bit = bitsPerm.at(i);
    //     LOG(INFO) << i << ") " << bit << " : " << histogram.bins().at(bit);
    // }

    ////////////////////////

    using Node = DecisionNode<BitSet, IndigoFingerprintTableView>;
    using Pair = std::pair<std::size_t, Node *>;
    
    Node &root = _decisionTree.getRoot();
    root.setInfo(IndigoFingerprintTableView(&_fingerprintTable));

    std::queue<Pair> nodes;
    nodes.push({0, &root});

    for (size_t nodeCount = 0; !nodes.empty(); nodeCount++) {

        Pair pair = nodes.front();

        std::size_t bit = pair.first;
        Node *node = pair.second;

        IndigoFingerprintTableView &view = node->getInfo();

        LOG(INFO) << "Bit: " << bit << ", size: " << view.size();

        if (bit < CHAR_BIT*qtr::IndigoFingerprint::sizeInBytes && view.size() > _maxLeafSize) {

            BitSet pred(_splittingStrategy(bit, view));
            typename Node::Children children = node->setPred(pred);

            auto viewPair = view.split2(pred);
            view = IndigoFingerprintTableView();

            children._next->setInfo(std::move(viewPair.first));
            children._false->setInfo(std::move(viewPair.second));

            nodes.push({bit + 1, children._next});
            nodes.push({bit + 1, children._false});
        }
        
        nodes.pop();
    }

}

template<class SplittingStrategy>
std::vector<IndigoMolecule> DecisionTreeSearchEngine<SplittingStrategy>::findOverMolecules(const IndigoQueryMolecule &mol)
{
    std::vector<IndigoMolecule> result;

    indigoAromatize(mol.id());

    QtrIndigoFingerprint fingerprint(mol, "sub");
    int bitsCount = fingerprint.countBits();

    qtr::IndigoFingerprint fp;
    fp.setBytes(fingerprint.data());

    std::vector<const IndigoFingerprintTableView *> views = _decisionTree.search(fp);

    for(const IndigoFingerprintTableView *view : views) {
        for (IndigoFingerprintTableView::IndexType idx : *view) {
            
            qtr::IndigoFingerprint f = _fingerprintTable.at(idx);
            f &= fp;
            
            if (f.count() != bitsCount)
                continue;
            
            const IndigoMolecule &molecule = _molecules.at(idx);
            IndigoSubstructureMatcher matcher = _indigoSessionPtr->substructureMatcher(molecule);
        
            if (!matcher.match(mol))
                continue;

            result.push_back(molecule);
        }
    }

    return result;
}

std::size_t SplittingStrategyOptimal::operator()(std::size_t bitIndex, const IndigoFingerprintTableView &view)
{
    Histogram histogram(CHAR_BIT*qtr::IndigoFingerprint::sizeInBytes);
    
    for(IndigoFingerprintTableView::IndexType index : view) {
        const qtr::IndigoFingerprint &fp = view.table()->at(index);
        for(size_t bit = 0; bit < fp.size(); bit++)
            histogram.add(bit, Histogram::CounterType(fp.test(bit)));
    }

    std::size_t result = std::size_t(-1);
    Histogram::CounterType half = Histogram::CounterType(view.size() / 2);
    Histogram::CounterType deviation = half;

    for(std::size_t bin = 0; bin < histogram.bins().size(); bin++) {
        
        Histogram::CounterType num = histogram.bins().at(bin);
        Histogram::CounterType dev = (num > half ? num - half : half - num);
        
        if (dev < deviation) {
            deviation = dev;
            result = bin;
        }
    }

    return result;
}

template class DecisionTreeSearchEngine<SplittingStrategyTrivial>;
template class DecisionTreeSearchEngine<SplittingStrategyOptimal>;

} // namespace qtr