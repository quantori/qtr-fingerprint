#include "ExperimentArgs.h"

#include <absl/flags/flag.h>
#include <absl/flags/parse.h>
#include <glog/logging.h>

#include <map>

#include "SearchEngineFactory.h"

namespace {

    const std::map<std::string, SearchEngineType> stringToSearchEngineType = {
            {"QtrRDKit",  SearchEngineType::QtrRDKit},
            {"QtrIndigo", SearchEngineType::QtrIndigo},
            {"RDKit",     SearchEngineType::RDKit},
            {"Indigo",    SearchEngineType::Indigo},
    };

    const std::string allSearchEngineTypes = [] {
        std::string res;
        for (auto it = stringToSearchEngineType.begin(); it != stringToSearchEngineType.end(); it++) {
            if (it != stringToSearchEngineType.begin()) {
                res += ", ";
            }
            res += it->first;
        }
        return res;
    }();

}
ABSL_FLAG(std::string, SearchEngineType, "",
          "Specifies the type of search engine to be tested. Possible options include: " + allSearchEngineTypes);

ABSL_FLAG(int, MaxResults, -1,
          "Sets the maximum number of results to be retrieved for each query. "
          "Specify a positive integer to limit the number of results.");

ABSL_FLAG(double, TimeLimit, -1,
          "Sets the time limit in seconds for each query. "
          "Specify a positive value to limit the duration of the query.");

ABSL_FLAG(std::string, QueriesFile, "",
          "Specifies the file containing the queries to be tested in the experiment.");

ABSL_FLAG(std::string, DatasetDir, "",
          "Specifies the directory containing the dataset. The dataset directory should consist of .csv files.");

ABSL_FLAG(std::string, StatisticsFile, "",
          "Specifies the file where experiment statistics will be written.");


ExperimentArgs::ExperimentArgs(int argc, char **argv) {
    absl::ParseCommandLine(argc, argv);

    auto searchEngineTypeIt = stringToSearchEngineType.find(absl::GetFlag(FLAGS_SearchEngineType));
    if (searchEngineTypeIt == stringToSearchEngineType.end()) {
        LOG(ERROR) << "Invalid value for SearchEngineType. Possible options: " << allSearchEngineTypes;
        exit(1);
    }
    searchEngineType = searchEngineTypeIt->second;

    maxResults = absl::GetFlag(FLAGS_MaxResults);
    if (maxResults <= 0) {
        LOG(ERROR) << "MaxResults must be a positive integer.";
        exit(1);
    }

    timeLimit = absl::GetFlag(FLAGS_TimeLimit);
    if (timeLimit <= 0) {
        LOG(ERROR) << "TimeLimit must be a positive value.";
        exit(1);
    }

    queriesFile = absl::GetFlag(FLAGS_QueriesFile);
    if (!std::filesystem::exists(queriesFile)) {
        LOG(ERROR) << "Queries file does not exist: " << queriesFile;
        exit(1);
    }

    datasetDir = absl::GetFlag(FLAGS_DatasetDir);
    if (!std::filesystem::exists(datasetDir)) {
        LOG(ERROR) << "Dataset directory does not exist: " << datasetDir;
        exit(1);
    }
    if (!std::filesystem::is_directory(datasetDir)) {
        LOG(ERROR) << "Specified dataset directory is not a directory: " << datasetDir;
        exit(1);
    }

    statisticsFile = absl::GetFlag(FLAGS_StatisticsFile);
    if (statisticsFile.empty()) {
        LOG(ERROR) << "Statistics file path cannot be empty.";
        exit(1);
    }
}
