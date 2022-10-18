#include "Utils.h"

#include "crow.h"
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"

#include <glog/logging.h>

using namespace qtr;

ABSL_FLAG(std::string, database_path, "",
          "Path to molecular database");

ABSL_FLAG(size_t, port, 8080,
          "Web service port");

int main(int argc, char *argv[]) {
    
    google::InitGoogleLogging(argv[0]);
    absl::ParseCommandLine(argc, argv);

    std::string databasePath = absl::GetFlag(FLAGS_database_path);
    auto servicePort = static_cast<size_t>(absl::GetFlag(FLAGS_port));

    emptyArgument(databasePath, "Please specify database_path option");

    crow::SimpleApp app;

    CROW_ROUTE(app, "/")([&databasePath](){
        return databasePath;
    });

    app.port(servicePort).multithreaded().run();

    return 0;
}