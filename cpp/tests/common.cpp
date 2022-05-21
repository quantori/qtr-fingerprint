#include "common.h"

#include "IndigoInChI.h"
#include "IndigoQueryMolecule.h"
#include "indigo.h"
#include "indigo-inchi.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <vector>

using namespace indigo_cpp;

void testSearchEngine(SearchEnginePtr searchEngine, IndigoSessionPtr indigoSession)
{
    const std::string databaseName = "119697";
    if (std::filesystem::exists(databaseName))
        searchEngine->build(databaseName);
    else
        searchEngine->build(databaseName + ".sdf");

    const std::string queriesFile = "119697_queries.txt";
    std::ifstream fin(queriesFile);

    IndigoInChI indigoInChi(indigoSession);

    std::string line;
    while (std::getline(fin, line))
    {
        std::istringstream iss(line);

        std::string smiles;
        iss >> smiles;

        std::vector<std::string> inchiKeys;

        while(true) {
            std::string inchiKey;
            iss >> inchiKey;
            if (!iss) break;
            inchiKeys.push_back(inchiKey);
        }

        int mol = indigoLoadQueryMoleculeFromString(smiles.c_str());
        auto queryMolecule = IndigoQueryMolecule(mol, indigoSession);
   
        std::vector<IndigoMolecule> result = searchEngine->findOverMolecules(queryMolecule);

        std::vector<std::string> inchiKeysResult;
        for(const IndigoMolecule &molecule : result) {
            //std::string smilesResult = molecule.canonicalSmiles();
            std::string inchi = indigoInChi.getInChI(molecule);
            std::string inchiKey = indigoInchiGetInchiKey(inchi.c_str());
            inchiKeysResult.push_back(inchiKey);
        }

        std::sort(inchiKeysResult.begin(), inchiKeysResult.end());

        compareTwoVectors(inchiKeys, inchiKeysResult);
   }

}