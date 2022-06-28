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
    Histogram histogram = view.histogram();
    return histogram.findClosestBin(Histogram::CounterType(view.size() / 2));
}

template class DecisionTreeSearchEngine<SplittingStrategyTrivial>;
template class DecisionTreeSearchEngine<SplittingStrategyOptimal>;

} // namespace qtr