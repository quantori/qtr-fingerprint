#include "BingoSearchEngine.h"

using namespace indigo_cpp;

BingoSearchEngine::BingoSearchEngine(const IndigoSessionPtr &indigoSessionPtr) : _indigoSessionPtr(indigoSessionPtr) {}

BingoSearchEngine::~BingoSearchEngine() {
    if (_db != -1)
        bingoCloseDatabase(_db);
}

bool static endsWith(const std::string &a, const std::string &b) {
    return a.size() >= b.size() &&
           (0 == a.compare(a.size() - b.size(), b.size(), b));
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
}

std::vector<indigo_cpp::IndigoMolecule>
BingoSearchEngine::findOverMolecules(const indigo_cpp::IndigoQueryMolecule &mol) {
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
