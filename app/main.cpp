#include "IndigoMolecule.h"
#include "IndigoSession.h"

#include <iostream>

using namespace indigo_cpp;


int main(void) {

    IndigoSessionPtr indigoSessionPtr = IndigoSession::create();
    IndigoMolecule molecule = indigoSessionPtr->loadMoleculeFromFile("query.mol");
    
    molecule.aromatize();
    std::cout << molecule.canonicalSmiles() << std::endl;

    return 0;
}