#include <glog/logging.h>

#include "BuildArgs.h"
#include "DatabaseBuilderFactory.h"
#include "Profiling.h"

using namespace std;
using namespace qtr;

int main(int argc, char *argv[]) {
    initLogging(argv, google::INFO, "build_db.info", true);
    BuildArgs args(argc, argv);

    try {
        auto builder = DatabaseBuilderFactory::create(args.dbType());
        builder->build(args);

        for (auto &[label, time]: ProfilingPool::getStatistics()) {
            LOG(INFO) << label << ": " << time;
        }
    }
    catch (std::exception &e) {
        LOG_ERROR_AND_EXIT(e.what());
    }

    return 0;
}