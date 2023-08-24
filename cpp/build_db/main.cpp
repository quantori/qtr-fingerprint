#include <glog/logging.h>

#include "BuildArgs.h"
#include "DatabaseBuilderFactory.h"

using namespace std;
using namespace qtr;

int main(int argc, char *argv[]) {
    initLogging(argv, google::INFO, "build_db.info", true);
    BuildArgs args(argc, argv);

    try {
        TimeMeasurer statisticCollector;
        auto builder = DatabaseBuilderFactory::create(args.dbType());
        builder->build(args, statisticCollector);

        for (auto &[label, time]: statisticCollector) {
            LOG(INFO) << label << ": " << time;
        }
    }
    catch (std::exception& e) {
        logErrorAndExit(e.what());
    }

    return 0;
}