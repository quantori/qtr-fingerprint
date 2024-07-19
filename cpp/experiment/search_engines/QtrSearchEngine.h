#pragma once

#include "BallTree.h"
#include "SearchEngineConcept.h"

template<SearchEngine SE>
class QtrSearchEngine {
public:
    using FingerprintType = SE::FingerprintType;
    using MoleculeType = SE::MoleculeType;
    using StorageMoleculeType = SE::StorageMoleculeType;
    using QueryMoleculeType = SE::QueryMoleculeType;


    explicit QtrSearchEngine(std::vector<std::string> &&smilesDataset) {
        std::vector<std::pair<std::unique_ptr<StorageMoleculeType>, std::unique_ptr<FingerprintType>>> data;
        for (auto &smiles: smilesDataset) {
            std::unique_ptr<MoleculeType> mol;
            try {
                mol = smilesToMolecule(smiles);
            } catch (const std::exception &e) {
                LOG(WARNING) << "Skip smiles: " << smiles << " error: " << e.what();
                continue;
            }
            if (mol == nullptr) {
                LOG(WARNING) << "Can't parse smiles: " << smiles;
                continue;
            }
            auto fp = std::make_unique<FingerprintType>(*mol);
            auto storageMol = moleculeToStorageMolecule(*mol);
            data.emplace_back(std::move(storageMol), std::move(fp));
        }
        _ballTree = std::make_unique<BallTreeType>(std::move(data));
    }

    explicit QtrSearchEngine(std::vector<std::unique_ptr<StorageMoleculeType>> &&moleculesDataset) {
        std::vector<std::pair<std::unique_ptr<StorageMoleculeType>, std::unique_ptr<FingerprintType>>> data;
        for (auto &storageMol: moleculesDataset) {
            auto mol = storageMoleculeToMolecule(*storageMol);
            auto fp = std::make_unique<FingerprintType>(*mol);
            data.emplace_back(std::move(storageMol), std::move(fp));
        }
        _ballTree = std::make_unique<BallTreeType>(std::move(data));
    }

    explicit QtrSearchEngine(
            std::vector<std::pair<std::unique_ptr<StorageMoleculeType>, std::unique_ptr<FingerprintType>>> &&data) {
        _ballTree = std::make_unique<BallTreeType>(std::move(data));
    }

    std::vector<uint64_t> getMatches(const QueryMoleculeType &queryMol, int maxResults, bool &stopFlag) {
        return _ballTree->getMatches(queryMol, maxResults, stopFlag);
    }

    std::vector<uint64_t>
    getMatches(const QueryMoleculeType &mol, const FingerprintType &fingerprint, int maxResults, bool &stopFlag) {
        return _ballTree->getMatches(mol, fingerprint, maxResults, stopFlag);
    }

    static std::unique_ptr<MoleculeType> smilesToMolecule(const std::string &smiles) {
        return SE::smilesToMolecule(smiles);
    }

    static std::unique_ptr<QueryMoleculeType> smilesToQueryMolecule(const std::string &smiles) {
        return SE::smilesToQueryMolecule(smiles);
    }

    static std::unique_ptr<MoleculeType> storageMoleculeToMolecule(const StorageMoleculeType &storageMolecule) {
        return SE::storageMoleculeToMolecule(storageMolecule);
    }

    static std::unique_ptr<StorageMoleculeType> moleculeToStorageMolecule(const MoleculeType &molecule) {
        return SE::moleculeToStorageMolecule(molecule);
    }

private:
    using BallTreeType = BallTree<SE>;

    std::unique_ptr<BallTreeType> _ballTree;
};
