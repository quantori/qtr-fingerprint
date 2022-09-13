#include "SearchEngineFactory.h"

#include "BingoSearchEngine.h"
#include "ExhaustiveSearchEngine.h"

namespace qtr {

    SearchEnginePtr SearchEngineFactory::create(
            SearchEngineType searchEngineType,
            indigo_cpp::IndigoSessionPtr indigoSessionPtr) {

        if (searchEngineType == SearchEngineType::BINGO)
            return std::make_shared<BingoSearchEngine>(indigoSessionPtr);

        if (searchEngineType == SearchEngineType::EXHAUSTIVE)
            return std::make_shared<ExhaustiveSearchEngine>(indigoSessionPtr);

        throw std::invalid_argument("Unknown search engine type!");
    }

} // namespace qtr
