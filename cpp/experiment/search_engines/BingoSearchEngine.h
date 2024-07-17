#pragma once

#include "IndigoFingerprint.h"

#include "BingoNoSQL.h"
#include "IndigoSession.h"
#include "indigo.h"

#include <filesystem>

class BingoSearchEngine {
public:
    using FingerprintType = IndigoFingerprint;
    using MoleculeType = indigo_cpp::IndigoMolecule;
    using QueryMoleculeType = indigo_cpp::IndigoQueryMolecule;

    explicit BingoSearchEngine(const std::vector<std::string> &smiles);

    explicit BingoSearchEngine(std::vector<std::unique_ptr<MoleculeType>> &&molecules);

    explicit BingoSearchEngine(
            std::vector<std::pair<std::unique_ptr<MoleculeType>, std::unique_ptr<FingerprintType>>> &&data);

    std::vector<uint64_t> getMatches(const QueryMoleculeType &queryMol, int maxResults, bool &stopFlag);

    std::vector<uint64_t>
    getMatches(const QueryMoleculeType &mol, const FingerprintType &fingerprint, int maxResults, bool &stopFlag);

    static std::unique_ptr<MoleculeType> smilesToMolecule(const std::string &smiles);

    static std::unique_ptr<QueryMoleculeType> smilesToQueryMolecule(const std::string &smiles);

    ~BingoSearchEngine();

private:

    BingoSearchEngine();

    std::filesystem::path _dbFilePath;
    indigo_cpp::BingoMolecule _db;
};
