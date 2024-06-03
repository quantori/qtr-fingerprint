#include <iostream>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/Fingerprints/Fingerprints.h>

int main(int argc, char **argv) {
  RDKit::ROMol *mol1 = RDKit::SmilesToMol("Cc1ccccc1");
  std::cout << "Number of atoms " << mol1->getNumAtoms() << std::endl;
  std::unique_ptr<ExplicitBitVect> mfp(RDKit::PatternFingerprintMol(*mol1));
}
