#include "CompleteSearchEngineProfiler.h"

using namespace indigo_cpp;

static long double percentHigh(std::vector<long double> durations, long double percent) {
    assert(!durations.empty());
    assert(percent >= 0 && percent <= 1);
    size_t count = ceil(durations.size() * percent);
    std::sort(durations.rbegin(), durations.rend());
    long double sum = 0;
    for (size_t i = 0; i < count; ++i) {
        sum += durations[i];
    }
    return sum / (long double)count;
}

void CompleteSearchEngineProfiler::profile(const std::string &path,
                                           SearchEngineInterface &searchEngine,
                                           const std::vector<IndigoQueryMolecule> &queries) {
    SearchEngineProfiler profiler(path, searchEngine);
    std::cout << std::setprecision(9) << std::fixed;

    std::cout << "Building search engine..." << std::endl;
    try {
        profiler.buildSearchEngine();
    } catch (...) {
        std::cerr << "Build failed!" << std::endl;
        throw;
    }
    auto buildDuration = profiler.getBuildDuration();
    std::cout << "Build completed" << std::endl;
    std::vector<long double> queryDurations(queries.size());
    if (!queries.empty()) {
        size_t maxW = std::to_string(queries.size() - 1).size();
        for (size_t i = 0; i < queries.size(); ++i) {
            try {
                profiler.profile(queries[i]);
            } catch (...) {
                std::cerr << "Query " << i << " failed!" << std::endl;
                throw;
            }
            std::cout << "Query " << std::setw(maxW) << i << " completed with time: "
                      << profiler.getLastProfileDuration().count() << std::endl;
            queryDurations[i] = profiler.getLastProfileDuration().count();
        }
    }
    auto durationWithoutBuild = profiler.getTotalDuration() - buildDuration;

    std::cout << "Profiling is completed" << std::endl;
    std::cout << "Build time:                " << buildDuration.count() << std::endl;
    if (!queries.empty()) {
        std::cout << "Average query time:        " << durationWithoutBuild.count() / double(queries.size())
                  << std::endl;
        std::cout << "10% high query time:       " << percentHigh(queryDurations, 0.1) << std::endl;
        std::cout << "1% high query time:        " << percentHigh(queryDurations, 0.01) << std::endl;
        std::cout << "0.1% high query time:      " << percentHigh(queryDurations, 0.001) << std::endl;
    }
    std::cout << "Total time(without build): " << durationWithoutBuild.count() << std::endl;
    std::cout << "Total time(with build):    " << profiler.getTotalDuration().count() << std::endl;
}