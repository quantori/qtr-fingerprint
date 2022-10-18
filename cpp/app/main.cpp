#include "Utils.h"

#include "crow.h"
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

#include <glog/logging.h>

using namespace qtr;

ABSL_FLAG(std::string, database_path, "",
          "Path to molecular database");

int main(int argc, char *argv[]) {
    
    google::InitGoogleLogging(argv[0]);
    absl::ParseCommandLine(argc, argv);

    std::string databasePath = absl::GetFlag(FLAGS_database_path);

    emptyArgument(databasePath, "Please specify database_path option");

    crow::SimpleApp app;

    CROW_ROUTE(app, "/")([&databasePath](){
        return databasePath;
    });

    app.port(18080).multithreaded().run();

    return 0;
}