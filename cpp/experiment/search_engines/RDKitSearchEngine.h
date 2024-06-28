#pragma once

#include "GraphMol/SubstructLibrary/SubstructLibrary.h"

#include "RDKitFingerprint.h"

#include <filesystem>

class RDKitSearchEngine {
public:
    using FingerprintType = RDKitFingerprint;

    explicit RDKitSearchEngine(const std::vector<std::string> &smiles);

    explicit RDKitSearchEngine(std::vector<std::unique_ptr<RDKit::ROMol>> &&molecules);

    explicit RDKitSearchEngine(
            std::vector<std::pair<std::unique_ptr<RDKit::ROMol>, std::unique_ptr<FingerprintType>>> &&data);

    std::vector<uint64_t> getMatches(const RDKit::ROMol &queryMol, int maxResults, bool &stopFlag);

    static std::unique_ptr<RDKit::ROMol> smilesToMolecule(const std::string &smiles);

private:
    std::shared_ptr<RDKit::SubstructLibrary> _substructLibrary;
};
