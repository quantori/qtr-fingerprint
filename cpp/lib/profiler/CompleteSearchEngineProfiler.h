#pragma once

#include "IndigoQueryMolecule.h"

#include "SearchEngineInterface.h"
#include "SearchEngineProfiler.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>

namespace CompleteSearchEngineProfiler {
    /**
     * Output's profiling to std::out
     *
     * @param path Path to file, on which the search engine will be build
     * @param searchEngine Implementation of SearchEngineInterface, not built yet
     * @param queries Vector of queries
     */
    void profile(const std::string &path,
                 SearchEngineInterface &searchEngine,
                 const std::vector<indigo_cpp::IndigoQueryMolecule> &queries);
};