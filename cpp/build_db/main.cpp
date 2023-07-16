#include <glog/logging.h>

#include "Args.h"
#include "build_types/buildQtr.h"
#include "build_types/buildBingoNoSQL.h"

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
    Args args(argc, argv);
    TimeMeasurer statisticCollector;

    if (args.dbType() == DatabaseType::QtrDrive ||
        args.dbType() == DatabaseType::QtrRam) {
        buildQtrDb(args, statisticCollector);
    } else if (args.dbType() == DatabaseType::BingoNoSQL) {
        buildBingoNoSQLDb(args, statisticCollector);
    }

    for (auto &[label, time]: statisticCollector) {
        LOG(INFO) << label << ": " << time;
    }

    return 0;
}