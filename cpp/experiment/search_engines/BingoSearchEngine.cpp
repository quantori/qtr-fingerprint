#include <glog/logging.h>
#include <execution>

#include "IndigoSearchEngine.h"
#include "indigo_molecule.h"
#include "indigo.h"
#include "bingo-nosql.h"
#include "GlobalIndigoSession.h"
#include "IndigoException.h"


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

IndigoSearchEngine::IndigoSearchEngine(std::unique_ptr<std::vector<std::string>> &&smiles) : IndigoSearchEngine() {
    std::for_each(std::execution::par, smiles->begin(), smiles->end(), [&](const std::string &smiles) {
        try {
            MoleculeType mol = globalIndigoSession->loadMolecule(smiles);
            mol.aromatize();
            {
                // no mutex as _db must be thread safe
                _db.insertRecord(mol);
            }
        }
        catch (const indigo_cpp::IndigoException &e) {
            LOG(ERROR) << "Indigo error: " << e.what();
        }
    });
}

IndigoSearchEngine::IndigoSearchEngine(
        std::unique_ptr<std::vector<std::unique_ptr<MoleculeType>>> &&molecules) : IndigoSearchEngine() {
    // no mutex as _db must be thread safe
    std::for_each(std::execution::par, molecules->begin(), molecules->end(), [&](std::unique_ptr<MoleculeType>& mol) {
        _db.insertRecord(*mol);
    });
}

IndigoSearchEngine::IndigoSearchEngine(
        std::unique_ptr<std::vector<std::pair<std::unique_ptr<MoleculeType>, std::unique_ptr<FingerprintType>>>> &&data)
        : IndigoSearchEngine() {
    // TODO: is it possible to do not ignore fingerprint?
    std::for_each(std::execution::par, data->begin(), data->end(), [&](decltype(data->operator[](0))& entity) {
        _db.insertRecord(*entity.first);
    });
}

IndigoSearchEngine::IndigoSearchEngine() : _dbFilePath(generateDBPath()), _db(
        indigo_cpp::BingoMolecule::createDatabaseFile(globalIndigoSession, _dbFilePath, "")) {
}

std::vector<uint64_t>
IndigoSearchEngine::getMatches(const IndigoSearchEngine::QueryMoleculeType &queryMol, int maxResults, bool &stopFlag) {
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

std::unique_ptr<IndigoSearchEngine::MoleculeType> IndigoSearchEngine::smilesToMolecule(const std::string &smiles) {
    auto mol = std::make_unique<MoleculeType>(globalIndigoSession->loadMolecule(smiles));
    mol->aromatize();
    return mol;
}

std::unique_ptr<IndigoSearchEngine::QueryMoleculeType>
IndigoSearchEngine::smilesToQueryMolecule(const std::string &smiles) {
    auto mol = std::make_unique<QueryMoleculeType>(globalIndigoSession->loadQueryMolecule(smiles));
    mol->aromatize();
    return mol;
}

std::vector<uint64_t> IndigoSearchEngine::getMatches(const IndigoSearchEngine::QueryMoleculeType &queryMol,
                                                     const IndigoSearchEngine::FingerprintType &fingerprint,
                                                     int maxResults, bool &stopFlag) {
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

IndigoSearchEngine::~IndigoSearchEngine() {
    std::filesystem::remove_all(_dbFilePath);
}

std::unique_ptr<IndigoSearchEngine::MoleculeType>
IndigoSearchEngine::storageMoleculeToMolecule(const IndigoSearchEngine::StorageMoleculeType &storageMolecule) {
    return std::make_unique<MoleculeType>(storageMolecule);
}

std::unique_ptr<IndigoSearchEngine::StorageMoleculeType>
IndigoSearchEngine::moleculeToStorageMolecule(const IndigoSearchEngine::MoleculeType &molecule) {
    return std::make_unique<StorageMoleculeType>(molecule);
}
