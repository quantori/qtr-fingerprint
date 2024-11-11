#pragma once

#include <execution>

#include "Trie.h"
#include "SearchEngineConcept.h"

template<SearchEngine SE>
class TrieSearchEngine {
public:
    using FingerprintType = typename SE::FingerprintType;
    using MoleculeType = typename SE::MoleculeType;
    using StorageMoleculeType = typename SE::StorageMoleculeType;
    using QueryMoleculeType = typename SE::QueryMoleculeType;

    explicit TrieSearchEngine(std::unique_ptr<std::vector<std::string>>&& smilesDataset, 
                            size_t maxDepth = std::numeric_limits<size_t>::max()) {
        auto data = std::make_unique<std::vector<std::pair<std::unique_ptr<StorageMoleculeType>, std::unique_ptr<FingerprintType>>>>();
        std::mutex dataMutex;
        std::for_each(std::execution::par, smilesDataset->begin(), smilesDataset->end(),
                      [&](const std::string& smiles) {
                          std::unique_ptr<MoleculeType> mol;
                          try {
                              mol = smilesToMolecule(smiles);
                          } catch (const std::exception& e) {
                              LOG(WARNING) << "Skip smiles: " << smiles << " error: " << e.what();
                              return;
                          }
                          if (mol == nullptr) {
                              LOG(WARNING) << "Can't parse smiles: " << smiles;
                              return;
                          }
                          auto fp = std::make_unique<FingerprintType>(*mol);
                          auto storageMol = moleculeToStorageMolecule(*mol);
                          {
                              std::lock_guard<std::mutex> lockGuard(dataMutex);
                              data->emplace_back(std::move(storageMol), std::move(fp));
                          }
                      });
        _trie = std::make_unique<TrieType>(std::move(data), maxDepth);
    }

    explicit TrieSearchEngine(std::unique_ptr<std::vector<std::unique_ptr<StorageMoleculeType>>>&& moleculesDataset,
                            size_t maxDepth = std::numeric_limits<size_t>::max()) {
        auto data = std::make_unique<std::vector<std::pair<std::unique_ptr<StorageMoleculeType>, std::unique_ptr<FingerprintType>>>>();
        std::mutex m;
        std::for_each(std::execution::par, moleculesDataset->begin(), moleculesDataset->end(),
                      [&](const auto& storageMol) {
                          auto mol = storageMoleculeToMolecule(*storageMol);
                          auto fp = std::make_unique<FingerprintType>(*mol);
                          {
                              std::lock_guard<std::mutex> lockGuard(m);
                              data->emplace_back(std::move(storageMol), std::move(fp));
                          }
                      });
        _trie = std::make_unique<TrieType>(std::move(data), maxDepth);
    }

    explicit TrieSearchEngine(
            std::unique_ptr<std::vector<std::pair<std::unique_ptr<StorageMoleculeType>, std::unique_ptr<FingerprintType>>>>&& data,
            size_t maxDepth = std::numeric_limits<size_t>::max()) {
        _trie = std::make_unique<TrieType>(std::move(data), maxDepth);
    }

    std::vector<uint64_t> getMatches(const QueryMoleculeType& queryMol, int maxResults, bool& stopFlag) {
        return _trie->getMatches(queryMol, maxResults, stopFlag);
    }

    std::vector<uint64_t>
    getMatches(const QueryMoleculeType& mol, const FingerprintType& fingerprint, int maxResults, bool& stopFlag) {
        return _trie->getMatches(mol, fingerprint, maxResults, stopFlag);
    }

    static std::unique_ptr<MoleculeType> smilesToMolecule(const std::string& smiles) {
        return SE::smilesToMolecule(smiles);
    }

    static std::unique_ptr<QueryMoleculeType> smilesToQueryMolecule(const std::string& smiles) {
        return SE::smilesToQueryMolecule(smiles);
    }

    static std::unique_ptr<MoleculeType> storageMoleculeToMolecule(const StorageMoleculeType& storageMolecule) {
        return SE::storageMoleculeToMolecule(storageMolecule);
    }

    static std::unique_ptr<StorageMoleculeType> moleculeToStorageMolecule(const MoleculeType& molecule) {
        return SE::moleculeToStorageMolecule(molecule);
    }

private:
    using TrieType = Trie<StorageMoleculeType, QueryMoleculeType, FingerprintType>;
    std::unique_ptr<TrieType> _trie;
}; 