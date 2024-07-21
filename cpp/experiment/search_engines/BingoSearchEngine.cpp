#include <glog/logging.h>

#include "BingoSearchEngine.h"
#include "indigo_molecule.h"
#include "indigo.h"
#include "bingo-nosql.h"
#include "GlobalIndigoSession.h"

#include "QtrBingoNoSQL.h"

namespace {
    std::mt19937 random_generator(0);

    std::filesystem::path generateDBPath() {
        std::filesystem::path res;
        do {
            uint32_t dirId = random_generator();
            std::string dirName = "bingoDb_" + std::to_string(dirId);
            res = std::filesystem::temp_directory_path() / dirName;
        } while (std::filesystem::exists(res));
//        LOG(INFO) << "Generated name: " << res;
        return res;
    }
}

BingoSearchEngine::BingoSearchEngine(std::unique_ptr<std::vector<std::string>> &&smiles) : BingoSearchEngine() {
    for (auto &s: *smiles) { // TODO: parallelize
        MoleculeType mol = globalIndigoSession->loadMolecule(s);
        _db.insertRecord(mol);
    }
}

BingoSearchEngine::BingoSearchEngine(
        std::unique_ptr<std::vector<std::unique_ptr<MoleculeType>>> &&molecules) : BingoSearchEngine() {
    for (auto &mol: *molecules) { // TODO: parallelize?
        _db.insertRecord(*mol);
    }
}

BingoSearchEngine::BingoSearchEngine(
        std::unique_ptr<std::vector<std::pair<std::unique_ptr<MoleculeType>, std::unique_ptr<FingerprintType>>>> &&data)
        : BingoSearchEngine() {
    // TODO: is it possible to do not ignore fingerprint?
    for (auto &[mol, _]: *data) { // TODO: parallelize?
        _db.insertRecord(*mol);
    }
}

BingoSearchEngine::BingoSearchEngine() : _dbFilePath(generateDBPath()), _db(
        QtrBingoNoSQL::createDatabaseFile(globalIndigoSession, _dbFilePath, "")) {

}

std::vector<uint64_t>
BingoSearchEngine::getMatches(const BingoSearchEngine::QueryMoleculeType &queryMol, int maxResults, bool &stopFlag) {
    std::vector<uint64_t> result;
    auto subMatcher = _db.searchSub(queryMol, "");
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

std::unique_ptr<BingoSearchEngine::QueryMoleculeType>
BingoSearchEngine::smilesToQueryMolecule(const std::string &smiles) {
    auto mol = std::make_unique<QueryMoleculeType>(globalIndigoSession->loadQueryMolecule(smiles));
    mol->aromatize();
    return mol;
}

std::vector<uint64_t> BingoSearchEngine::getMatches(const BingoSearchEngine::QueryMoleculeType &queryMol,
                                                    const BingoSearchEngine::FingerprintType &fingerprint,
                                                    int maxResults, bool &stopFlag) {
    std::vector<uint64_t> result;
    auto subMatcher = _db.searchSub(queryMol, fingerprint.array(), "");
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

BingoSearchEngine::~BingoSearchEngine() {
    std::filesystem::remove_all(_dbFilePath);
}

std::unique_ptr<BingoSearchEngine::MoleculeType>
BingoSearchEngine::storageMoleculeToMolecule(const BingoSearchEngine::StorageMoleculeType &storageMolecule) {
    return std::make_unique<MoleculeType>(storageMolecule);
}

std::unique_ptr<BingoSearchEngine::StorageMoleculeType>
BingoSearchEngine::moleculeToStorageMolecule(const BingoSearchEngine::MoleculeType &molecule) {
    return std::make_unique<StorageMoleculeType>(molecule);
}
