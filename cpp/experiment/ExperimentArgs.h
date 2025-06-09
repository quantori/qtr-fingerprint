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
    int ballTreeDepth; // -1 means calculate depth automatically
    double fingerprintRatio;

    ExperimentArgs(int argc, char *argv[]);
};
