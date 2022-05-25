#include "common.h"

#include "IndigoInChI.h"
#include "IndigoQueryMolecule.h"
#include "indigo.h"
#include "indigo-inchi.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using namespace indigo_cpp;
using namespace std;

static pair<string, vector<string>> parseLine(const string &line)
{
    istringstream iss(line);

    string smiles;
    iss >> smiles;

    vector<string> inchiKeys;

    while(true) {
        string inchiKey;
        iss >> inchiKey;
        if (!iss) break;
        inchiKeys.push_back(inchiKey);
    }

    return {smiles, inchiKeys};
}

static vector<string> parseResult(const vector<IndigoMolecule> &result, IndigoInChI &indigoInChi)
{
    vector<string> inchiKeys;

    for(const IndigoMolecule &molecule : result) {
        std::string inchi = indigoInChi.getInChI(molecule);
        std::string inchiKey = indigoInchiGetInchiKey(inchi.c_str());
        inchiKeys.push_back(inchiKey);
    }

    sort(inchiKeys.begin(), inchiKeys.end());

    return inchiKeys;
}

void testSearchEngine(
    SearchEnginePtr searchEngine,
    IndigoSessionPtr indigoSession,
    const std::string &fileSdf,
    const std::string &fileQueries)
{
    IndigoInChI indigoInChi(indigoSession);

    searchEngine->build(fileSdf);

    std::ifstream fin(fileQueries);

    std::string line;
    while (std::getline(fin, line))
    {
        auto [smiles, inchiKeys] = parseLine(line);

        int mol = indigoLoadQueryMoleculeFromString(smiles.c_str());
        auto queryMolecule = IndigoQueryMolecule(mol, indigoSession);
   
        std::vector<IndigoMolecule> result = searchEngine->findOverMolecules(queryMolecule);

        auto inchiKeysResult = parseResult(result, indigoInChi);

        compareTwoVectors(inchiKeys, inchiKeysResult);
   }    
}

std::filesystem::path getDataDir()
{
    using namespace std::filesystem;
    const path currentDir = testing::UnitTest::GetInstance()->original_working_dir();
    const path dataDir = currentDir / path("./data");
    return dataDir;
}
