#include "SearchEngineProfiler.h"

using namespace indigo_cpp;

SearchEngineProfiler::SearchEngineProfiler(const std::string &path, SearchEngineInterface &searchEngine)
        : path(path), searchEngine(searchEngine) {}


void SearchEngineProfiler::buildSearchEngine() {
    auto startTime = std::chrono::high_resolution_clock::now();

    this->searchEngine.build(this->path);

    auto endTime = std::chrono::high_resolution_clock::now();
    lastDuration = endTime - startTime;
}

std::vector<IndigoMolecule> SearchEngineProfiler::profile(const IndigoQueryMolecule &query) {
    auto startTime = std::chrono::high_resolution_clock::now();

    auto result = searchEngine.findOverMolecules(query);

    auto endTime = std::chrono::high_resolution_clock::now();
    lastDuration = endTime - startTime;
    totalDuration += lastDuration;

    return result;
}

std::chrono::duration<long double> SearchEngineProfiler::getLastDuration() const noexcept {
    return lastDuration;
}

std::chrono::duration<long double> SearchEngineProfiler::getTotalDuration() const noexcept {
    return totalDuration;
}