#include "RunArgs.h"

using namespace std;

ABSL_FLAG(string, dbName, "",
          "Name of folders with data base's files");

ABSL_FLAG(string, dbType, "",
          "Possible types: "
                  FLAG_NAME(QtrDrive) ", "
                  FLAG_NAME(QtrRam) ", "
                  FLAG_NAME(BingoNoSQL));

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

ABSL_FLAG(std::string, summaryFile, "",
          "File where summary of the run should be stored (for " FLAG_NAME(FromFile) " mode only)");

