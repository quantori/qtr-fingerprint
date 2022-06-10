#include "SearchEngineFactory.h"

#include "BingoSearchEngine.h"
#include "ExhaustiveSearchEngine.h"

#include <exception>

qtr::SearchEnginePtr qtr::SearchEngineFactory::create(SearchEngineType searchEngineType, indigo_cpp::IndigoSessionPtr indigoSessionPtr) {
    if (searchEngineType == SearchEngineType::BINGO)    
        return std::make_shared<BingoSearchEngine>(indigoSessionPtr);
    else if (searchEngineType == SearchEngineType::EXHAUSTIVE)
        return std::make_shared<ExhaustiveSearchEngine>(indigoSessionPtr);
    else {
        throw std::invalid_argument("Unknown search engine type!");
        return nullptr;
    }
}
