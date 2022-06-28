#include "SearchEngineFactory.h"

#include "BingoSearchEngine.h"
#include "DecisionTreeSearchEngine.h"
#include "ExhaustiveSearchEngine.h"

#include <exception>

namespace qtr {

SearchEnginePtr SearchEngineFactory::create(
    SearchEngineType searchEngineType,
    indigo_cpp::IndigoSessionPtr indigoSessionPtr,
    size_t decisionTreeMaxLeafSize) {

    if (searchEngineType == SearchEngineType::BINGO)    
        return std::make_shared<BingoSearchEngine>(indigoSessionPtr);
    
    if (searchEngineType == SearchEngineType::EXHAUSTIVE)
        return std::make_shared<ExhaustiveSearchEngine>(indigoSessionPtr);
    
    if (searchEngineType == SearchEngineType::DECISION_TREE_TRIVIAL)
        return std::make_shared<DecisionTreeSearchEngine<SplittingStrategyTrivial>>(indigoSessionPtr, decisionTreeMaxLeafSize);

    if (searchEngineType == SearchEngineType::DECISION_TREE_OPTIMAL)
        return std::make_shared<DecisionTreeSearchEngine<SplittingStrategyOptimal>>(indigoSessionPtr, decisionTreeMaxLeafSize);
    
    throw std::invalid_argument("Unknown search engine type!");
    return nullptr;
}

} // namespace qtr
