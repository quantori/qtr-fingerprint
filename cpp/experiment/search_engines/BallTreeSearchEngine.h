#pragma once

#include <execution>

#include "BallTree.h"
#include "SearchEngineConcept.h"

template<SearchEngine SE>
class BallTreeSearchEngine {
public:
    using FingerprintType = SE::FingerprintType;
    using MoleculeType = SE::MoleculeType;
    using StorageMoleculeType = SE::StorageMoleculeType;
    using QueryMoleculeType = SE::QueryMoleculeType;


    explicit BallTreeSearchEngine(std::unique_ptr<std::vector<std::string>> &&smilesDataset) {
        auto data = std::make_unique<std::vector<std::pair<std::unique_ptr<StorageMoleculeType>, std::unique_ptr<FingerprintType>>>>();
        std::mutex dataMutex;
        std::for_each(std::execution::par, smilesDataset->begin(), smilesDataset->end(),
                      [&](const std::string &smiles) {
                          std::unique_ptr<MoleculeType> mol;
                          try {
                              mol = smilesToMolecule(smiles);
                          } catch (const std::exception &e) {
                              LOG(WARNING) << "Skip smiles: " << smiles << " error: " << e.what();
                              return;
                          }
                          if (mol == nullptr) {
                              LOG(WARNING) << "Can't parse smiles: " << smiles;
                              return;
                          }
                          auto fp = std::make_unique<FingerprintType>(*mol);
                          assert(fp != nullptr);
                          auto storageMol = moleculeToStorageMolecule(*mol);
                          {
                              std::lock_guard<std::mutex> lockGuard(dataMutex);
                              data->emplace_back(std::move(storageMol), std::move(fp));
                          }
                      });
        _ballTree = std::make_unique<BallTreeType>(std::move(data));
    }

    explicit BallTreeSearchEngine(
            std::unique_ptr<std::vector<std::unique_ptr<StorageMoleculeType>>> &&moleculesDataset) {
        auto data = std::make_unique<std::vector<std::pair<std::unique_ptr<StorageMoleculeType>, std::unique_ptr<FingerprintType>>>>;
        std::mutex m;
        std::for_each(std::execution::par, moleculesDataset->begin(), moleculesDataset->end(),
                      [&](const StorageMoleculeType &storageMol) {
                          auto mol = storageMoleculeToMolecule(*storageMol);
                          auto fp = std::make_unique<FingerprintType>(*mol);
                          {
                              std::lock_guard<std::mutex> lockGuard(m);
                              data->emplace_back(std::move(storageMol), std::move(fp));
                          }
                      });
        _ballTree = std::make_unique<BallTreeType>(std::move(data));
    }

    explicit BallTreeSearchEngine(
            std::unique_ptr<std::vector<std::pair<std::unique_ptr<StorageMoleculeType>, std::unique_ptr<FingerprintType>>>> &&data) {
        _ballTree = std::make_unique<BallTreeType>(std::move(data));
    }

    std::vector<uint64_t> getMatches(QueryMoleculeType &queryMol, int maxResults, bool &stopFlag) {
        return _ballTree->getMatches(queryMol, maxResults, stopFlag);
    }

    std::vector<uint64_t>
    getMatches(QueryMoleculeType &mol, const FingerprintType &fingerprint, int maxResults, bool &stopFlag) {
        return _ballTree->getMatches(mol, fingerprint, maxResults, stopFlag);
    }

    static std::unique_ptr<MoleculeType> smilesToMolecule(const std::string &smiles) {
        return SE::smilesToMolecule(smiles);
    }

    static std::unique_ptr<QueryMoleculeType> smilesToQueryMolecule(const std::string &smiles) {
        return SE::smilesToQueryMolecule(smiles);
    }

    static std::unique_ptr<MoleculeType> storageMoleculeToMolecule(StorageMoleculeType &storageMolecule) {
        return SE::storageMoleculeToMolecule(storageMolecule);
    }

    static std::unique_ptr<StorageMoleculeType> moleculeToStorageMolecule(MoleculeType &molecule) {
        return SE::moleculeToStorageMolecule(molecule);
    }

private:
    using BallTreeType = BallTree<SE>;

    std::unique_ptr<BallTreeType> _ballTree;
};
