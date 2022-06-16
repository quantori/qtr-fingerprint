#include "SearchEngineFactory.h"

#include "BingoSearchEngine.h"
#include "DecisionTreeSearchEngine.h"
#include "ExhaustiveSearchEngine.h"

#include <exception>

namespace qtr {

SearchEnginePtr SearchEngineFactory::create(SearchEngineType searchEngineType, indigo_cpp::IndigoSessionPtr indigoSessionPtr) {
    if (searchEngineType == SearchEngineType::BINGO)    
        return std::make_shared<BingoSearchEngine>(indigoSessionPtr);
    else if (searchEngineType == SearchEngineType::EXHAUSTIVE)
        return std::make_shared<ExhaustiveSearchEngine>(indigoSessionPtr);
    else if (searchEngineType == SearchEngineType::DECISION_TREE)
        return std::make_shared<DecisionTreeSearchEngine>(indigoSessionPtr);
    else {
        throw std::invalid_argument("Unknown search engine type!");
        return nullptr;
    }
}

} // namespace qtr
