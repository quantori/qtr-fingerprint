#include "SearchEngineFactory.h"

#include <stdexcept>

SearchEngine SearchEngineFactory::create(SearchEngineType type) {
    throw std::logic_error("SearchEngine::create is not implemented");
}
