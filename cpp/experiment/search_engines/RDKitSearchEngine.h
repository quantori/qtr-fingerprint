#pragma once

#include "GraphMol/SubstructLibrary/SubstructLibrary.h"

#include "RDKitFingerprint.h"

#include <filesystem>

class RDKitSearchEngine {
public:
    using FingerprintType = RDKitFingerprint;
    using MoleculeType = RDKit::ROMol;
    using QueryMoleculeType = RDKit::ROMol;

    explicit RDKitSearchEngine(const std::vector<std::string> &smiles);

    explicit RDKitSearchEngine(std::vector<std::unique_ptr<MoleculeType>> &&molecules);

    explicit RDKitSearchEngine(
            std::vector<std::pair<std::unique_ptr<MoleculeType >, std::unique_ptr<FingerprintType>>> &&data);

    std::vector<uint64_t> getMatches(const QueryMoleculeType &queryMol, int maxResults, bool &stopFlag);

    static std::unique_ptr<MoleculeType> smilesToMolecule(const std::string &smiles);

    static std::unique_ptr<MoleculeType> smilesToQueryMolecule(const std::string &smiles);

private:
    std::shared_ptr<RDKit::SubstructLibrary> _substructLibrary;
};
