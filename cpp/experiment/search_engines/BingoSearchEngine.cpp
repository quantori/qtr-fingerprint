#include "BingoSearchEngine.h"
#include "indigo_molecule.h"
#include "indigo.h"
#include "bingo-nosql.h"
#include "GlobalIndigoSession.h"

namespace {
    std::mt19937 random_generator(0);
}

BingoSearchEngine::BingoSearchEngine(const std::vector<std::string> &smiles) {
    for (auto &s: smiles) {
        MoleculeType mol = globalIndigoSession->loadMolecule(s);
        _db->insertRecord(mol);
    }
}

BingoSearchEngine::BingoSearchEngine(std::vector<std::unique_ptr<MoleculeType>> &&molecules) : BingoSearchEngine() {
    for (auto &mol: molecules) {
        _db->insertRecord(*mol);
    }
}

BingoSearchEngine::BingoSearchEngine(
        std::vector<std::pair<std::unique_ptr<MoleculeType>, std::unique_ptr<FingerprintType>>> &&data) {
    // TODO: is it possible to do not ignore fingerprint?
    for (auto &[mol, _]: data) {
        _db->insertRecord(*mol);
    }
}

BingoSearchEngine::BingoSearchEngine() {
    do {
        uint32_t dirId = random_generator();
        std::string dirName = "bingoDb_" + std::to_string(dirId);
        _dbDir = std::filesystem::temp_directory_path() / dirName;
    } while (std::filesystem::exists(_dbDir));
    _db = std::make_unique<indigo_cpp::BingoMolecule>(
            indigo_cpp::BingoMolecule::createDatabaseFile(globalIndigoSession, _dbDir, ""));
}

std::vector<uint64_t>
BingoSearchEngine::getMatches(const BingoSearchEngine::QueryMoleculeType &queryMol, int maxResults, bool &stopFlag) {
    // TODO: check molecule aromatized?
    std::vector<uint64_t > result;
    auto subMatcher = _db->searchSub(queryMol, "");
    for (auto &mol: subMatcher) {
        if (stopFlag) {
            break;
        }
        result.push_back(mol.getId()); // TODO: Fix abstract Ids?
        if (result.size() >= maxResults) {
            break;
        }
    }
    return result;
}

std::unique_ptr<BingoSearchEngine::MoleculeType> BingoSearchEngine::smilesToMolecule(const std::string &smiles) {
    auto mol = std::make_unique<MoleculeType>(globalIndigoSession->loadMolecule(smiles));
    mol->aromatize();
    return mol;
}
