#include "SearchEngineProfiler.h"

using namespace indigo_cpp;

qtr::SearchEngineProfiler::SearchEngineProfiler(const std::string &path, SearchEngineInterface &searchEngine)
        : _path(path), _searchEngine(searchEngine) {}


void qtr::SearchEngineProfiler::buildSearchEngine() {
    auto startTime = std::chrono::high_resolution_clock::now();

    _searchEngine.build(_path);

    auto endTime = std::chrono::high_resolution_clock::now();
    _buildDuration = endTime - startTime;
    _totalDuration += _buildDuration;
}

std::vector<IndigoMolecule> qtr::SearchEngineProfiler::profile(const IndigoQueryMolecule &query) {
    auto startTime = std::chrono::high_resolution_clock::now();

    auto result = _searchEngine.findOverMolecules(query);

    auto endTime = std::chrono::high_resolution_clock::now();
    _lastProfileDuration = endTime - startTime;
    _totalDuration += _lastProfileDuration;

    return result;
}

const std::chrono::duration<long double> &qtr::SearchEngineProfiler::getBuildDuration() const {
    return _buildDuration;
}

std::chrono::duration<long double> qtr::SearchEngineProfiler::getLastProfileDuration() const noexcept {
    return _lastProfileDuration;
}

std::chrono::duration<long double> qtr::SearchEngineProfiler::getTotalDuration() const noexcept {
    return _totalDuration;
}
