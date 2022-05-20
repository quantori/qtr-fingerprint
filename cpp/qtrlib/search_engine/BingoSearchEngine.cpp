#include "BingoSearchEngine.h"

using namespace indigo_cpp;

BingoSearchEngine::BingoSearchEngine(const IndigoSessionPtr &indigoSessionPtr) : _indigoSessionPtr(indigoSessionPtr) {}

BingoSearchEngine::~BingoSearchEngine() {
    if (currentDatabaseState != NOT_CREATED && currentDatabaseState != CLOSED) {
        bingoCloseDatabase(_db);
        currentDatabaseState = CLOSED;
    }
}

void BingoSearchEngine::build(const std::string &path) {
    if (endsWith(path, ".sdf")) {
        std::string dbName = path.substr(0, path.size() - strlen(".sdf"));
        _db = bingoCreateDatabaseFile(dbName.c_str(), "molecule", "");
        int iterator = indigoIterateSDFile(path.c_str());
        bingoInsertIteratorObj(_db, iterator);
    } else {
        _db = bingoLoadDatabaseFile(path.c_str(), "");
    }
    currentDatabaseState = LOADED;
}

std::vector<indigo_cpp::IndigoMolecule>
BingoSearchEngine::findOverMolecules(const indigo_cpp::IndigoQueryMolecule &mol) {
    if (currentDatabaseState != LOADED)
        throw std::runtime_error("Database is not loaded");
    std::vector<IndigoMolecule> result;
    int resultIterator = _indigoSessionPtr->_checkResult(bingoSearchSub(_db, mol.id(), ""));
    int currentId = bingoGetObject(resultIterator);
    do {
        IndigoMolecule curr = _indigoSessionPtr->loadMolecule(indigoCanonicalSmiles(currentId));
        result.emplace_back(curr);
    } while (bingoNext(resultIterator));
    bingoEndSearch(resultIterator);
    return result;
}
