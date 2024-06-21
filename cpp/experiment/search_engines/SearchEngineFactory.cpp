#include "SearchEngineFactory.h"

#include <stdexcept>

std::unique_ptr<SearchEngine> SearchEngineFactory::create(SearchEngineType type, const std::filesystem::path &datasetDir) {
    switch (type) {
        case SearchEngineType::QtrRDKit:
            throw std::invalid_argument("QtrRDKit search engine is not implemented");
        case SearchEngineType::QtrIndigo:
            throw std::invalid_argument("QtrIndigo search engine is not implemented");
        case SearchEngineType::RDKit:
            throw std::invalid_argument("RDKit search engine is not implemented");
        case SearchEngineType::Indigo:
            throw std::invalid_argument("Indigo search engine is not implemented");
        default:
            throw std::invalid_argument("Invalid SearchEngineType");
    }
}
