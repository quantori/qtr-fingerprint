#include "indigo.h"

#include "IndigoMolecule.h"
#include "IndigoQueryMolecule.h"
#include "IndigoSession.h"
#include "IndigoSDFileIterator.h"
#include "IndigoException.h"
#include "BingoNoSQL.h"
#include "BingoResultIterator.h"

#include <iostream>
#include <cstring>
#include <cassert>
#include <memory>
#include <vector>

using namespace indigo_cpp;
using namespace std;

/*
 * ./program_name DB_PATH load FILENAME.sdf
 *  to load database from sdf file
 *
 * ./program_name DB_PATH
 * to profile database
 */

int main(int argc, char **argv) {
    assert(argc > 1 && "No database path");
    string databasePath = argv[1];
    IndigoSessionPtr indigoSessionPtr = IndigoSession::create();
    if (argc > 2 && strcmp(argv[2], "load") == 0) {
        assert(argc > 3 && "No input file");
        auto db = BingoMolecule::createDatabaseFile(indigoSessionPtr, databasePath);
        db.insertIterator(indigoSessionPtr->iterateSDFile(argv[3]));
        cout << "Done" << endl;
        return 0;
    }
    try {
        auto db = BingoMolecule::loadDatabaseFile(indigoSessionPtr, databasePath);
        cout << "query:\n";
        auto queryMolBase = indigoSessionPtr->loadMoleculeFromFile("query.mol");
        auto query = indigoSessionPtr->loadQueryMolecule(queryMolBase.canonicalSmiles());
        cout << queryMolBase.canonicalSmiles() << '\n';
        cout << "result:\n";

        int i = 0;
        std::unique_ptr<IndigoMolecule> molPtr;

        for (BingoResult<IndigoMolecule> &u : db.searchSub(query)) {
            if (!molPtr) molPtr = make_unique<IndigoMolecule>(u.getTarget());
            cout << i++ << ") " << u.getId() << " " << molPtr->canonicalSmiles() << endl;
        }
    } catch(exception &e) {
        cout << "Error!\n" << e.what() << endl;
    }
    return 0;
}
