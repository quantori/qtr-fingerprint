/**
 * All flags are defined in one place to avoid multiple definition when one flag is used in multiple applications
 */

#include "absl/flags/flag.h"

#include <string>

#include "ArgsBase.h"

using namespace std;

ABSL_FLAG(string, dbType, "",
          "Possible types: "
                  FLAG_NAME(QtrDrive) ", "
                  FLAG_NAME(QtrRam) ", "
                  FLAG_NAME(BingoNoSQL));

ABSL_FLAG(string, dbName, "",
          "Name of folders with database's files");

ABSL_FLAG(bool, properties, true,
          "True if properties should be concerned. False otherwise");


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

// preprocessing args:
ABSL_FLAG(string, preprocessingType, "",
          "Source files type. "
          "Possible types: "
                  FLAG_NAME(SDF) ", "
                  FLAG_NAME(CSV));

ABSL_FLAG(string, destDir, "", "Destination directory where preprocessed files should be stored");

ABSL_FLAG(string, targetFilesType, "",
          "Target files type (for " FLAG_NAME(SDF) " preprocessing type only). "
                                                   "Possible types: "
                  FLAG_NAME(RawBucket) ", "
                  FLAG_NAME(Tables));


ABSL_FLAG(vector<string>, dataDirs, {},
          "Path to directories where data are stored");

ABSL_FLAG(string, otherDataDir, "",
          "Path to directory with other files");

ABSL_FLAG(uint64_t, threads, -1,
          "Number of threads to process leafs.");

ABSL_FLAG(string, mode, "",
          "Possible modes: "
                  FLAG_NAME(Interactive) ", "
                  FLAG_NAME(FromFile) ", "
                  FLAG_NAME(Web));

ABSL_FLAG(string, queriesFile, "",
          "File to load test molecules from it (for " FLAG_NAME(FromFile) " mode only)");

ABSL_FLAG(uint64_t, ansCount, -1,
          "The maximum number of answers to be found");

ABSL_FLAG(double, timeLimit, -1,
          "Single request time limit in seconds");

ABSL_FLAG(string, summaryFile, "",
          "File where summary of the run should be stored (for " FLAG_NAME(FromFile) " mode only)");




