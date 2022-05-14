#include "BingoSearchEngine.h"

using namespace indigo_cpp;

BingoSearchEngine::BingoSearchEngine() {
    indigoSetSessionId(sessionId);
}

BingoSearchEngine::~BingoSearchEngine() {
    if (_db != -1)
        bingoCloseDatabase(_db);
}

bool static endsWith(const std::string &a, const std::string &b) {
    return a.size() >= b.size() &&
           (0 == a.compare(a.size() - b.size(), b.size(), b));
}

#include <iostream> // TODO for debug

void BingoSearchEngine::build(const std::string &path) {
    indigoSetSessionId(sessionId);
    if (endsWith(path, ".sdf")) {
        std::string dbName = path.substr(0, path.size() - strlen(".sdf"));

        int db = bingoCreateDatabaseFile(dbName.c_str(), "molecule", "");
        int iterator = indigoIterateSDFile(path.c_str());
        bingoInsertIteratorObj(db, iterator);
    } else {
        _db = bingoLoadDatabaseFile(path.c_str(), "");
    }
}

std::vector<indigo_cpp::IndigoMolecule>
BingoSearchEngine::findOverMolecules(const indigo_cpp::IndigoQueryMolecule &mol) {
    indigoSetSessionId(sessionId);
//    IndigoSessionPtr indigoSessionPtr(new IndigoSession(sessionId)); // uncomment this line, and it wont work anymore
    std::vector<IndigoMolecule> result;
    int resultIterator = bingoSearchSub(_db, mol.id(), "");
    if (resultIterator == -1) {
        std::cerr << indigoGetLastError() << '\n';
        exit(-1);
    }
    int currentId = bingoGetObject(resultIterator);
    do {
//            auto curr = indigoSessionPtr->loadMolecule(indigoCanonicalSmiles(currentId));
//            std::cout << curr.canonicalSmiles() << '\n'; // ideally, this should work

        std::cout << indigoCanonicalSmiles(currentId) << '\n';
    } while (bingoNext(resultIterator));

    bingoEndSearch(resultIterator);
    return result;
}

