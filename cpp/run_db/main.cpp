#include <string>

#include <glog/logging.h>

#include "Utils.h"

#include "RunArgs.h"
#include "SearchDataLoader.h"
#include "RunModeFactory.h"

using namespace std;
using namespace qtr;

int main(int argc, char *argv[]) {
    initLogging(argv, google::INFO, "run_db.info", true);
    RunArgs args(argc, argv);

    try {
        auto searchData = SearchDataLoader::load(args);
        unique_ptr<RunMode> mode = RunModeFactory::create(args, searchData);
        mode->run();
    } catch (const exception &e) {
        LOG_ERROR_AND_EXIT(e.what());
    }
    return 0;
}