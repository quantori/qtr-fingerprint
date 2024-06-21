#include <glog/logging.h>
#include "ExperimentArgs.h"


int main(int argc, char *argv[]) {
    // Init logging
    google::InitGoogleLogging(argv[0]);
    google::SetLogDestination(google::INFO, "experiment.info");
    FLAGS_alsologtostderr = true;

    ExperimentArgs args(argc, argv);



    return 0;
}