#include "indigo.h"

#include "IndigoMolecule.h"
#include "IndigoSession.h"
#include "IndigoSDFileIterator.h"

#include <fstream>
#include <iostream>

using namespace indigo_cpp;
using namespace std;

int main()
{
    IndigoSessionPtr indigoSessionPtr = IndigoSession::create();
    IndigoSDFileIterator fileIt = indigoSessionPtr->iterateSDFile("../test1/119697.sdf");

    ofstream fout("../test1/fingerprints_119697.txt");

    for(IndigoMoleculeSPtr mol : fileIt) {
        const char *id = indigoGetProperty(mol->id(), "PUBCHEM_COMPOUND_CID");
        fout << id << " ";
        mol->aromatize();
        int fingerprint = indigoFingerprint(mol->id(), "sub");
        const char *fp = indigoToString(fingerprint);
        fout << fp << endl;
    }

    return 0;
}