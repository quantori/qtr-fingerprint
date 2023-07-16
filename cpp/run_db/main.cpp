#include <string>

#include <glog/logging.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>

#include "Utils.h"

#include "Args.h"
#include "TimeTicker.h"
#include "search_data/SearchDataLoader.h"
#include "modes/RunMode.h"
#include "modes/InteractiveMode.h"
#include "modes/web/WebMode.h"
#include "modes/FromFileMode.h"

using namespace std;
using namespace qtr;

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
          "File to load test molecules from");

ABSL_FLAG(uint64_t, ansCount, -1,
          "The maximum number of answers to be found");

ABSL_FLAG(double, timeLimit, -1,
          "Single request time limit in seconds");

int main(int argc, char *argv[]) {
    initLogging(argv, google::INFO, "run_db.info", true);
    Args args(argc, argv);

    try {
        TimeTicker timeTicker;
        auto searchData = SearchDataLoader::load(args, timeTicker);
        timeTicker.tick("Db data loading");

        shared_ptr<RunMode> mode = nullptr;
        if (args.mode() == Args::Mode::Interactive)
            mode = make_shared<InteractiveMode>(searchData);
        else if (args.mode() == Args::Mode::FromFile)
            mode = make_shared<FromFileMode>(searchData, args.queriesFile());
        else if (args.mode() == Args::Mode::Web)
            mode = make_shared<WebMode>(searchData);
        mode->run();
    } catch (const exception &e) {
        logErrorAndExit(e.what());
    }
    return 0;
}