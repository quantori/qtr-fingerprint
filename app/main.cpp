#include "IndigoException.h"
#include "IndigoMolecule.h"
#include "IndigoSession.h"

#include <iostream>

using namespace indigo_cpp;


int main(void) {

    IndigoSessionPtr indigoSessionPtr = IndigoSession::create();
    try {
        IndigoMolecule molecule = indigoSessionPtr->loadMoleculeFromFile("query.mol");
        std::cout << molecule.canonicalSmiles() << std::endl;
    } catch(IndigoException &ie) {
        std::cout << ie.what() << std::endl;
    }

    return 0;
}