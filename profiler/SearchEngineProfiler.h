#pragma once

#include "IndigoMolecule.h"
#include "IndigoQueryMolecule.h"

#include "SearchEngineInterface.h"

#include <chrono>
#include <vector>

class SearchEngineProfiler {
public:
    /**
     * Build search engine profiler from path and implementation of SearchEngineInterface
     *
     * @param path Path to the file with .sdf.gz extension
     */
    SearchEngineProfiler(const std::string &path, SearchEngineInterface &searchEngine);

    ~SearchEngineProfiler() = default;

    /**
     * Build search engine with the path given in constructor
     *
     * Duration of build also calculated, but wont be counted in total time
     */
    void buildSearchEngine();

    /**
     * Profile search engine on query
     *
     * @param query Query indigo molecule
     * @return Vector of parent molecules, to get time of operation call getLastDuration
     */
    std::vector<indigo_cpp::IndigoMolecule> profile(const indigo_cpp::IndigoQueryMolecule &query);

    /**
     * @return last profile duration
     */
    [[nodiscard]] std::chrono::duration<long double> getLastDuration() const noexcept;

    /**
     * @return total profile duration
     */
    [[nodiscard]] std::chrono::duration<long double> getTotalDuration() const noexcept;

private:

    std::string path;
    SearchEngineInterface &searchEngine;
    std::chrono::duration<long double> lastDuration;
    std::chrono::duration<long double> totalDuration = std::chrono::duration<long double>::zero();
};