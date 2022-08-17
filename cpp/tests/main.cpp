#include <glog/logging.h>
#include <gtest/gtest.h>

#include <absl/flags/flag.h>
#include <absl/flags/parse.h>

#include <stdexcept>

#include "utils/GTestEnvironment.h"

ABSL_FLAG(std::string, data_dir_path, "",
          "Path to data dir");

ABSL_FLAG(std::string, big_data_dir_path, "",
          "Path to big data dir");

ABSL_FLAG(std::string, tmp_data_dir_path, "",
          "Path to temporal data dir");

int main(int argc, char **argv)
{
    google::InitGoogleLogging(argv[0]);
    ::testing::InitGoogleTest(&argc, argv);

    absl::ParseCommandLine(argc, argv);
    std::string dataDirPath = absl::GetFlag(FLAGS_data_dir_path);
    std::string bigDataDirPath = absl::GetFlag(FLAGS_big_data_dir_path);
    std::string tmpDataDirPath = absl::GetFlag(FLAGS_tmp_data_dir_path);

    ::testing::AddGlobalTestEnvironment(new qtr::GTestEnvironment(dataDirPath, bigDataDirPath, tmpDataDirPath));
    return RUN_ALL_TESTS();
}
