#include "SearchEngineProfiler.h"

using namespace indigo_cpp;

SearchEngineProfiler::SearchEngineProfiler(const std::string &path, SearchEngineInterface &searchEngine)
        : _path(path), _searchEngine(searchEngine) {}


void SearchEngineProfiler::buildSearchEngine() {
    auto startTime = std::chrono::high_resolution_clock::now();

    _searchEngine.build(_path);

    auto endTime = std::chrono::high_resolution_clock::now();
    _lastDuration = endTime - startTime;
}

std::vector<IndigoMolecule> SearchEngineProfiler::profile(const IndigoQueryMolecule &query) {
    auto startTime = std::chrono::high_resolution_clock::now();

    auto result = _searchEngine.findOverMolecules(query);

    auto endTime = std::chrono::high_resolution_clock::now();
    _lastDuration = endTime - startTime;
    _totalDuration += _lastDuration;

    return result;
}

std::chrono::duration<long double> SearchEngineProfiler::getLastDuration() const noexcept {
    return _lastDuration;
}

std::chrono::duration<long double> SearchEngineProfiler::getTotalDuration() const noexcept {
    return _totalDuration;
}