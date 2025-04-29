#pragma once

#include <string>
#include <filesystem>
#include "search/engines/SearchEngineInterface.h"
#include "search/engines/SearchEngineType.h"

struct ExperimentArgs {
    SearchEngineType searchEngineType;
    int maxResults;
    double timeLimit;
    std::filesystem::path queriesFile;
    std::filesystem::path datasetDir;
    std::filesystem::path queriesStatisticFile;
    std::filesystem::path searchEngineStatisticFile;

    ExperimentArgs(int argc, char *argv[]);
};
