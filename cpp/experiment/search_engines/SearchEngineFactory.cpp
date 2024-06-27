#include "SearchEngineFactory.h"

#include <stdexcept>

#include "RDKitSearchEngine.h"

std::unique_ptr<SearchEngine>
SearchEngineFactory::create(SearchEngineType type, const std::filesystem::path &datasetDir) {
    switch (type) {
        case SearchEngineType::QtrRDKit:
            throw std::invalid_argument("QtrRDKit search engine is not implemented");
        case SearchEngineType::QtrIndigo:
            throw std::invalid_argument("QtrIndigo search engine is not implemented");
        case SearchEngineType::RDKit:
            return std::make_unique<RDKitSearchEngine>(datasetDir);
        case SearchEngineType::Indigo:
            throw std::invalid_argument("Indigo search engine is not implemented");
        default:
            throw std::invalid_argument("Invalid SearchEngineType");
    }
}
