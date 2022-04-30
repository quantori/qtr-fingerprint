#pragma once

#include "IndigoMolecule.h"

#include "SearchEngineInterface.h"
#include "SearchEngineProfiler.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <vector>

namespace CompleteSearchEngineProfiler {
    /**
     * Output's profiling to std::out
     *
     * @param path Path to file, on which the search engine will be build
     * @param searchEngine Implementation of SearchEngineInterface, not built yet
     * @param queries Vector of queries
     */
    void profile(std::string path,
                        SearchEngineInterface &searchEngine,
                        std::vector<indigo_cpp::IndigoMolecule> queries);
};