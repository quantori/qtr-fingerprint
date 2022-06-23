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
    : FingerprintTableSearchEngine(indigoSessionPtr)
    , _maxLeafSize(maxLeafSize)
{}

template<class SplittingStrategy>
void DecisionTreeSearchEngine<SplittingStrategy>::build(const std::string &path)
{
    FingerprintTableSearchEngine::build(path);

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

            auto viewMap = view.split<bool, qtr::IndigoFingerprint>(pred);
            view = IndigoFingerprintTableView();

            children._next->setInfo(std::move(viewMap[true]));
            children._false->setInfo(std::move(viewMap[false]));

            nodes.push({bit + 1, children._next});
            nodes.push({bit + 1, children._false});
        }
        
        nodes.pop();
    }
}

template<class SplittingStrategy>
std::vector<const IndigoFingerprintTableView *> 
DecisionTreeSearchEngine<SplittingStrategy>::findTableViews(const qtr::IndigoFingerprint &fp) const
{
    return _decisionTree.search(fp);
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

        if (bitIndex == 0) {
            LOG(INFO) << "[" << bin << "] : " << num;
        }
        
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