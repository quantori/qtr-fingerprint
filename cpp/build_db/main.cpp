#include <glog/logging.h>

#include "BuildArgs.h"
#include "database_builders/DatabaseBuilderFactory.h"

using namespace std;
using namespace qtr;

ABSL_FLAG(string, dbType, "",
          "Possible types: "
                  FLAG_NAME(QtrDrive) ", "
                  FLAG_NAME(QtrRam) ", "
                  FLAG_NAME(BingoNoSQL));

ABSL_FLAG(string, dbName, "",
          "Name of folders with database's files");


ABSL_FLAG(string, sourceDir, "",
          "Path to directory where source data are stored");

ABSL_FLAG(vector<string>, destDirs, {},
          "Path to directories where data should be stored");


ABSL_FLAG(string, otherDestDir, "",
          "Path to directory where other data should be stored");

ABSL_FLAG(uint64_t, parallelizeDepth, 0,
          "Depth on which subtree parallelization starts");

ABSL_FLAG(uint64_t, treeDepth, 0,
          "ball tree depth");


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