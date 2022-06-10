#include "SearchEngineProfiler.h"

using namespace indigo_cpp;

namespace qtr {

SearchEngineProfiler::SearchEngineProfiler(const std::string &path, SearchEngineInterface &searchEngine)
        : _path(path), _searchEngine(searchEngine) {}


void SearchEngineProfiler::buildSearchEngine() {
    auto startTime = std::chrono::high_resolution_clock::now();

    _searchEngine.build(_path);

    auto endTime = std::chrono::high_resolution_clock::now();
    _buildDuration = endTime - startTime;
    _totalDuration += _buildDuration;
}

std::vector<IndigoMolecule> SearchEngineProfiler::profile(const IndigoQueryMolecule &query) {
    auto startTime = std::chrono::high_resolution_clock::now();

    auto result = _searchEngine.findOverMolecules(query);

    auto endTime = std::chrono::high_resolution_clock::now();
    _lastProfileDuration = endTime - startTime;
    _totalDuration += _lastProfileDuration;

    return result;
}

const std::chrono::duration<long double> &SearchEngineProfiler::getBuildDuration() const {
    return _buildDuration;
}

std::chrono::duration<long double> SearchEngineProfiler::getLastProfileDuration() const noexcept {
    return _lastProfileDuration;
}

std::chrono::duration<long double> SearchEngineProfiler::getTotalDuration() const noexcept {
    return _totalDuration;
}

} // namespace qtr
