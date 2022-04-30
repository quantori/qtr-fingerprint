#include "CompleteSearchEngineProfiler.h"

using namespace indigo_cpp;

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
    auto buildDuration = profiler.getLastDuration();
    std::cout << "Build completed" << std::endl;
    std::cout << "Build time: " << buildDuration.count() << std::endl;
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
                      << profiler.getLastDuration().count() << std::endl;
        }
    }

    std::cout << "Profiling is completed" << std::endl;
    std::cout << "Total time(without build): " << profiler.getTotalDuration().count() << std::endl;
    std::cout << "Total time(with build):    " << (profiler.getTotalDuration() + buildDuration).count() << std::endl;
}