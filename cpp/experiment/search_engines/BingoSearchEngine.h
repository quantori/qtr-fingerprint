#pragma once

#include "IndigoFingerprint.h"

#include "BingoNoSQL.h"
#include "IndigoSession.h"
#include "indigo.h"
//#include "../../third_party/Indigo/api/cpp/src/IndigoQueryMolecule.h"

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

    static std::unique_ptr<MoleculeType> smilesToMolecule(const std::string &smiles);

private:

    BingoSearchEngine();

    std::unique_ptr<indigo_cpp::BingoMolecule> _db;
    std::filesystem::path _dbDir;
};
