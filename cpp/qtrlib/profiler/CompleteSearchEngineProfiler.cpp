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
    return sum / (long double) count;
}

void qtr::CompleteSearchEngineProfiler::profile(const std::string &path,
                                           SearchEngineInterface &searchEngine,
                                           const std::vector<IndigoQueryMolecule> &queries) {
    LOG(INFO) << "Start profiling " << queries.size() << " molecules on search engine with path=" << path
              << std::setprecision(9) << std::fixed;
    SearchEngineProfiler profiler(path, searchEngine);

    try {
        profiler.buildSearchEngine();
    } catch (const std::exception &e) {
        LOG(ERROR) << "Build search engine is failed with message: " << e.what();
        throw;
    }
    auto buildDuration = profiler.getBuildDuration();
    std::vector<long double> queryDurations(queries.size());
    if (!queries.empty()) {
        size_t maxW = std::to_string(queries.size() - 1).size();
        for (size_t i = 0; i < queries.size(); ++i) {
            try {
                profiler.profile(queries[i]);
            } catch (const std::exception &e) {
                LOG(ERROR) << "Query #" << i << " failed with message: " << e.what();
                throw;
            }
            LOG(INFO) << "Query " << std::setw(maxW) << i << " completed with time: "
                      << profiler.getLastProfileDuration().count() << "s";
            queryDurations[i] = profiler.getLastProfileDuration().count();
        }
    }
    auto durationWithoutBuild = profiler.getTotalDuration() - buildDuration;

    LOG(INFO) << "Profiling is completed";
    LOG(INFO) << "Build time:                " << buildDuration.count() << "s";
    if (!queries.empty()) {
        LOG(INFO) << "Average query time:        " << durationWithoutBuild.count() / double(queries.size())
                  << "s";
        LOG(INFO) << "10% high query time:       " << percentHigh(queryDurations, 0.1) << "s";
        LOG(INFO) << "1% high query time:        " << percentHigh(queryDurations, 0.01) << "s";
        LOG(INFO) << "0.1% high query time:      " << percentHigh(queryDurations, 0.001) << "s";
    }
    LOG(INFO) << "Total time(without build): " << durationWithoutBuild.count() << "s";
    LOG(INFO) << "Total time(with build):    " << profiler.getTotalDuration().count() << "s";
}