#pragma once

#include "IndigoMolecule.h"
#include "IndigoQueryMolecule.h"

#include "SearchEngineInterface.h"

#include <chrono>
#include <vector>

namespace qtr {

class SearchEngineProfiler {
public:
    /**
     * Build search engine profiler from path and implementation of SearchEngineInterface
     *
     * @param path Path to the file with .sdf.gz extension
     * @param searchEngine SearchEngineInterface implementation
     */
    SearchEngineProfiler(const std::string &path, SearchEngineInterface &searchEngine);

    ~SearchEngineProfiler() = default;

    /**
     * Build search engine with the path given in constructor
     *
     * @see getBuildDuration
     */
    void buildSearchEngine();

    /**
     * Profile search engine on query
     *
     * @param query Query indigo molecule
     * @return Vector of parent molecules
     * @see getLastProfileDuration
     */
    std::vector<indigo_cpp::IndigoMolecule> profile(const indigo_cpp::IndigoQueryMolecule &query);

    /**
     * @return build engine duration
     */
    [[nodiscard]] const std::chrono::duration<long double> &getBuildDuration() const;

    /**
     * @return last profile duration
     */
    [[nodiscard]] std::chrono::duration<long double> getLastProfileDuration() const noexcept;

    /**
     * Total duration is a sum of build duration and all profile queries duration
     *
     * @return total profile duration
     */
    [[nodiscard]] std::chrono::duration<long double> getTotalDuration() const noexcept;

private:

    std::string _path;
    SearchEngineInterface &_searchEngine;
    std::chrono::duration<long double> _buildDuration;
    std::chrono::duration<long double> _lastProfileDuration;
    std::chrono::duration<long double> _totalDuration = std::chrono::duration<long double>::zero();
};

} // namespace qtr