#pragma once

#include "BallTree.h"
#include "RDKitSearchEngine.h"
#include "RDKitFingerprint.h"

class QtrRDKitSearchEngine {
public:
    explicit QtrRDKitSearchEngine(const std::vector<std::string> &smilesDataset);

    explicit QtrRDKitSearchEngine(std::vector<std::unique_ptr<RDKit::ROMol>> &&moleculesDataset);

    explicit QtrRDKitSearchEngine(
            std::vector<std::pair<std::unique_ptr<RDKit::ROMol>, std::unique_ptr<RDKitFingerprint>>> &&data);

    std::vector<uint64_t> getMatches(const RDKit::ROMol &queryMol, int maxResults, bool &stopFlag);

    static std::unique_ptr<RDKit::ROMol> smilesToMolecule(const std::string &smiles);

    using FingerprintType = RDKitSearchEngine::FingerprintType;

private:
    using BallTreeType = BallTree<RDKit::ROMol, RDKitFingerprint, RDKitSearchEngine>;

    std::unique_ptr<BallTreeType> _ballTree;
};
