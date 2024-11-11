#pragma once

#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <stdexcept>
#include <optional>
#include <iostream>

#include "FingerprintConcept.h"
#include "IndigoQueryMolecule.h"
#include "IndigoMolecule.h"
#include "IndigoFingerprint.h"
#include "IndigoSearchEngine.h"
#include "GlobalIndigoSession.h"
#include "IndigoSubstructureMatcher.h"
#include "IndigoUtils.h"

class IndigoBruteForceSearchEngine {
public:
    using FingerprintType = IndigoFingerprint;
    using MoleculeType = indigo_cpp::IndigoMolecule;
    using StorageMoleculeType = indigo_cpp::IndigoMolecule;
    using QueryMoleculeType = indigo_cpp::IndigoQueryMolecule;

    explicit IndigoBruteForceSearchEngine(std::unique_ptr<std::vector<std::string>> &&smiles) {
        for (const auto &smile: *smiles) {
            auto molecule = smilesToMolecule(smile);
            auto storageMolecule = moleculeToStorageMolecule(*molecule);
            auto fingerprint = std::make_unique<FingerprintType>(*molecule);
            molecules.emplace_back(std::move(storageMolecule));
            fingerprints.emplace_back(std::move(fingerprint));
        }
    }

    explicit IndigoBruteForceSearchEngine(std::unique_ptr<std::vector<std::unique_ptr<StorageMoleculeType>>> &&mols) {
        for (auto &storageMol: *mols) {
            auto mol = storageMoleculeToMolecule(*storageMol);
            auto fingerprint = std::make_unique<FingerprintType>(*mol);
            molecules.emplace_back(std::move(mol));
            fingerprints.emplace_back(std::move(fingerprint));
        }
    }

    explicit IndigoBruteForceSearchEngine(
            std::unique_ptr<std::vector<std::pair<std::unique_ptr<StorageMoleculeType>, std::unique_ptr<FingerprintType>>>> &&data
    ) {
        for (auto &[mol, fp]: *data) {
            molecules.emplace_back(std::move(mol));
            fingerprints.emplace_back(std::move(fp));
        }
    }

    std::vector<uint64_t> getMatches(const QueryMoleculeType &queryMol, int maxResults, bool &stopFlag) {
        FingerprintType queryFingerprint(queryMol);
        return getMatches(queryMol, queryFingerprint, maxResults, stopFlag);
    }

    std::vector<uint64_t>
    getMatches(const QueryMoleculeType &queryMol, const FingerprintType &queryFingerprint, int maxResults,
               bool &stopFlag) {
        std::vector<uint64_t> matches;
        for (size_t i = 0; i < fingerprints.size() && matches.size() < (size_t) maxResults && !stopFlag; i++) {
            if (queryFingerprint.isSubFingerprintOf(*fingerprints[i]) && isSubstructure(queryMol, *molecules[i])) {
                matches.push_back(i);
            }
        }
        return matches;
    }

    static std::unique_ptr<MoleculeType> smilesToMolecule(const std::string &smiles) {
        return IndigoSearchEngine::smilesToMolecule(smiles);
    }

    static std::unique_ptr<QueryMoleculeType> smilesToQueryMolecule(const std::string &smiles) {
        return IndigoSearchEngine::smilesToQueryMolecule(smiles);
    }

    static std::unique_ptr<MoleculeType> storageMoleculeToMolecule(const StorageMoleculeType &storageMolecule) {
        return IndigoSearchEngine::storageMoleculeToMolecule(storageMolecule);
    }

    static std::unique_ptr<StorageMoleculeType> moleculeToStorageMolecule(const MoleculeType &molecule) {
        return IndigoSearchEngine::moleculeToStorageMolecule(molecule);
    }

private:
    std::vector<std::unique_ptr<StorageMoleculeType>> molecules;
    std::vector<std::unique_ptr<FingerprintType>> fingerprints;
};
