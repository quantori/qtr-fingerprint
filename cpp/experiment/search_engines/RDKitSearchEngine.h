#pragma once

#include "SearchEngine.h"
#include "GraphMol/SubstructLibrary/SubstructLibrary.h"

#include <filesystem>

class RDKitSearchEngine : public SearchEngine {
public:
    RDKitSearchEngine(const std::filesystem::path &datasetDir);

    std::vector<uint64_t> getMatches(const std::string &querySmiles, int maxResults, bool &stopFlag) override;

private:
    std::shared_ptr<RDKit::SubstructLibrary> _substructLibrary;
};
