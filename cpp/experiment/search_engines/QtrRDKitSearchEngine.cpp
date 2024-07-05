#include "QtrRDKitSearchEngine.h"

std::vector<uint64_t> QtrRDKitSearchEngine::getMatches(const RDKit::ROMol &queryMol, int maxResults, bool &stopFlag) {
    return _ballTree->getMatches(queryMol, maxResults, stopFlag);
}

std::unique_ptr<RDKit::ROMol> QtrRDKitSearchEngine::smilesToMolecule(const std::string &smiles) {
    return RDKitSearchEngine::smilesToMolecule(smiles);
}

QtrRDKitSearchEngine::QtrRDKitSearchEngine(const std::vector<std::string> &smilesDataset) {
    std::vector<std::pair<std::unique_ptr<RDKit::ROMol>, std::unique_ptr<RDKitFingerprint>>> data;
    for (auto &smiles: smilesDataset) {
        std::unique_ptr<RDKit::ROMol> mol;
        try {
            mol = smilesToMolecule(smiles);
        } catch (const std::exception& e) {
            LOG(WARNING) << "Skip smiles: " << smiles << " error: " << e.what();
            continue;
        }
        if (mol == nullptr) {
            LOG(WARNING) << "Can't parse smiles: " << smiles;
            continue;
        }
        auto fp = std::make_unique<RDKitFingerprint>(*mol);
        data.emplace_back(std::move(mol), std::move(fp));
    }
    _ballTree = std::make_unique<BallTreeType>(std::move(data));
}

QtrRDKitSearchEngine::QtrRDKitSearchEngine(std::vector<std::unique_ptr<RDKit::ROMol>> &&moleculesDataset) {
    std::vector<std::pair<std::unique_ptr<RDKit::ROMol>, std::unique_ptr<RDKitFingerprint>>> data;
    for (auto &mol: moleculesDataset) {
        auto fp = std::make_unique<RDKitFingerprint>(*mol);
        data.emplace_back(std::move(mol), std::move(fp));
    }
    _ballTree = std::make_unique<BallTreeType>(std::move(data));
}

QtrRDKitSearchEngine::QtrRDKitSearchEngine(
        std::vector<std::pair<std::unique_ptr<RDKit::ROMol>, std::unique_ptr<RDKitFingerprint>>> &&data) {
    _ballTree = std::make_unique<BallTreeType>(std::move(data));
}
