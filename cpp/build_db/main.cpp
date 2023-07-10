#include <glog/logging.h>

#include "Args.h"
#include "buildQtr.h"
#include "buildBingoNoSQL.h"

using namespace std;
using namespace qtr;

ABSL_FLAG(std::string, dbType, "",
          "Possible types: "
                  FLAG_NAME(QtrDrive) ", "
                  FLAG_NAME(QtrRam) ", "
                  FLAG_NAME(BingoNoSQL));

ABSL_FLAG(std::string, dbName, "",
          "Name of folders with data base's files");


ABSL_FLAG(std::string, sourceDirPath, "",
          "Path to directory where source data are stored");

ABSL_FLAG(std::vector<std::string>, destDirPaths, {},
          "Path to directories where data should be stored");


ABSL_FLAG(std::string, otherDestDirPath, "",
          "Path to directory where other data should be stored");

ABSL_FLAG(uint64_t, parallelizeDepth, 0,
          "Depth on which subtree parallelization starts");

ABSL_FLAG(uint64_t, treeDepth, 0,
          "ball tree depth");


int main(int argc, char *argv[]) {
    initLogging(argv, google::INFO, "build_db.info", true);
    Args args(argc, argv);
    TimeMeasurer statisticCollector;

    if (args.dbType() == ArgsBase::DataBaseType::QtrDrive ||
        args.dbType() == ArgsBase::DataBaseType::QtrRam) {
        buildQtrDb(args, statisticCollector);
    } else if (args.dbType() == ArgsBase::DataBaseType::BingoNoSQL) {
        buildBingoNoSQLDb(args, statisticCollector);
    }

    for (auto &[label, time]: statisticCollector) {
        LOG(INFO) << label << ": " << time;
    }

    return 0;
}