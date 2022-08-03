#include <glog/logging.h>
#include <gtest/gtest.h>

#include <absl/flags/flag.h>
#include <absl/flags/parse.h>

#include "utils/GTestEnvironment.h"

ABSL_FLAG(std::string, data_dir_path, "",
          "Path to data dir");

ABSL_FLAG(std::string, big_data_dir_path, "",
          "Path to big data dir");

int main(int argc, char **argv)
{
    google::InitGoogleLogging(argv[0]);
    ::testing::InitGoogleTest(&argc, argv);

    absl::ParseCommandLine(argc, argv);
    std::string bigDataDirPath = absl::GetFlag(FLAGS_big_data_dir_path);
    std::string dataDirPath = absl::GetFlag(FLAGS_data_dir_path);

    ::testing::AddGlobalTestEnvironment(new qtr::GTestEnvironment(dataDirPath, bigDataDirPath));
    return RUN_ALL_TESTS();
}
