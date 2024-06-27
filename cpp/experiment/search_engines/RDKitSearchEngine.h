#pragma once

#include "GraphMol/SubstructLibrary/SubstructLibrary.h"

#include <filesystem>

class RDKitSearchEngine {
public:
    RDKitSearchEngine(const std::filesystem::path &datasetDir);

    std::vector<uint64_t> getMatches(const std::string &querySmiles, int maxResults, bool &stopFlag);

private:
    std::shared_ptr<RDKit::SubstructLibrary> _substructLibrary;
};
