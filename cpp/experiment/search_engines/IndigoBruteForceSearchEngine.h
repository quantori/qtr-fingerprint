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
    using MoleculeType = indigo::Molecule;
    using StorageMoleculeType = indigo::Molecule;
    using QueryMoleculeType = indigo::QueryMolecule;

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

    std::vector<uint64_t> getMatches(QueryMoleculeType &queryMol, int maxResults, bool &stopFlag) {
        FingerprintType queryFingerprint(queryMol);
        return getMatches(queryMol, queryFingerprint, maxResults, stopFlag);
    }

    std::vector<uint64_t>
    getMatches(QueryMoleculeType &queryMol, const FingerprintType &queryFingerprint, int maxResults,
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
        BufferScanner bufferScanner(smiles.c_str(), smiles.size(), false);
        SmilesLoader loader(bufferScanner);
        auto mol = std::make_unique<MoleculeType>();
        loader.loadMolecule(*mol);
        return mol;
    }

    static std::unique_ptr<QueryMoleculeType> smilesToQueryMolecule(const std::string &smiles) {
        BufferScanner bufferScanner(smiles.c_str(), smiles.size(), false);
        SmilesLoader loader(bufferScanner);
        auto mol = std::make_unique<QueryMoleculeType>();
        loader.loadQueryMolecule(*mol);
        return mol;
    }

    static std::unique_ptr<MoleculeType> storageMoleculeToMolecule(StorageMoleculeType &storageMolecule) {
        auto mol = std::make_unique<MoleculeType>();
        storageMolecule.clone(*mol);
        return mol;
    }

    static std::unique_ptr<StorageMoleculeType> moleculeToStorageMolecule(MoleculeType &molecule) {
        auto storageMol = std::make_unique<MoleculeType>();
        molecule.clone(*storageMol);
        return storageMol;
    }

private:
    std::vector<std::unique_ptr<StorageMoleculeType>> molecules;
    std::vector<std::unique_ptr<FingerprintType>> fingerprints;
};
