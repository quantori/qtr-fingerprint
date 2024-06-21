#include <glog/logging.h>
#include <absl/flags/flag.h>
#include <absl/flags/parse.h>

#include <map>

static std::map<std::string,

ABSL_FLAG(std::string, SearchEngineType, "", "Search engine type to be tested");



int main(int argc, char * argv[]) {
    // Init logging
    google::InitGoogleLogging(argv[0]);
    google::SetLogDestination(google::INFO, "experiment.info");
    FLAGS_alsologtostderr = true;

    // Parse command line arguments
    absl::ParseCommandLine(argc, argv);


    return 0;
}