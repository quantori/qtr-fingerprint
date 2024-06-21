#pragma once

#include <string>
#include <filesystem>

#include "SearchEngineFactory.h"

struct ExperimentArgs {
    SearchEngineType searchEngineType;
    int maxResults;
    double timeLimit;
    std::filesystem::path queriesFile;
    std::filesystem::path datasetDir;
    std::filesystem::path statisticsFile;

    ExperimentArgs(int argc, char *argv[]);
};
