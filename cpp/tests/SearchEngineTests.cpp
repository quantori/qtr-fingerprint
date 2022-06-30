#include "SearchEngineTests.h"
#include "Common.h"
#include "utils/DataPathManager.h"

#include "IndigoInChI.h"
#include "IndigoQueryMolecule.h"
#include "indigo.h"
#include "indigo-inchi.h"

#include <glog/logging.h>

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using namespace indigo_cpp;
using namespace qtr;
using namespace std;
using namespace std::filesystem;

static pair<string, vector<string>> parseQueryLine(const string &line)
{
    istringstream iss(line);

    string smiles;
    iss >> smiles;

    vector<string> inchiKeys;

    while (true)
    {
        string inchiKey;
        iss >> inchiKey;
        if (!iss)
            break;
        inchiKeys.push_back(inchiKey);
    }

    return {smiles, inchiKeys};
}

static vector<string> parseResult(const vector<IndigoMolecule> &result, IndigoInChI &indigoInChi)
{
    vector<string> inchiKeys;

    for (const IndigoMolecule &molecule : result)
    {
        std::string inchi = indigoInChi.getInChI(molecule);
        std::string inchiKey = indigoInchiGetInchiKey(inchi.c_str());
        inchiKeys.push_back(inchiKey);
    }

    sort(inchiKeys.begin(), inchiKeys.end());

    return inchiKeys;
}

static void testSearchEngine(
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
        auto [smiles, inchiKeys] = parseQueryLine(line);
        LOG(INFO) << smiles;

        int mol = indigoLoadQueryMoleculeFromString(smiles.c_str());
        auto queryMolecule = IndigoQueryMolecule(mol, indigoSession);

        std::vector<IndigoMolecule> result = searchEngine->findOverMolecules(queryMolecule);

        auto inchiKeysResult = parseResult(result, indigoInChi);

        compareTwoVectors(inchiKeys, inchiKeysResult);
    }
}

static void tesBuildSearchEngine(
    SearchEnginePtr searchEngine,
    const std::string &fileSdf)
{
    searchEngine->build(fileSdf);
}

namespace qtr
{

    void SearchEngineTests::testPubchem10()
    {
        testSearchEngine(
            _searchEnginePtr,
            _indigoSessionPtr,
            _dataDir / path("pubchem_10.sdf"),
            _dataDir / path("pubchem_10_queries.txt"));
    }

    void SearchEngineTests::testPubchem100()
    {
        testSearchEngine(
            _searchEnginePtr,
            _indigoSessionPtr,
            _dataDir / path("pubchem_100.sdf"),
            _dataDir / path("pubchem_100_queries.txt"));
    }

    void SearchEngineTests::testPubchem994()
    {
        testSearchEngine(
            _searchEnginePtr,
            _indigoSessionPtr,
            _dataDir / path("pubchem_994.sdf"),
            _dataDir / path("pubchem_994_queries.txt"));
    }

    void SearchEngineTests::testPubchem119697()
    {
        testSearchEngine(
            _searchEnginePtr,
            _indigoSessionPtr,
            DataPathManager::getBigDataDir() / path("pubchem_119697.sdf"),
            DataPathManager::getBigDataDir() / path("pubchem_119697_queries.txt"));
    }

    void SearchEngineTests::testBuildPubchem119697()
    {
        tesBuildSearchEngine(
            _searchEnginePtr,
            DataPathManager::getBigDataDir() / path("pubchem_119697.sdf"));
    }

    void SearchEngineTests::testPubchem300000()
    {
        testSearchEngine(
            _searchEnginePtr,
            _indigoSessionPtr,
            DataPathManager::getBigDataDir() / path("pubchem_300000.sdf"),
            DataPathManager::getBigDataDir() / path("pubchem_300000_queries.txt"));
    }

    void SearchEngineTests::testBuildPubchem300000()
    {
        tesBuildSearchEngine(
            _searchEnginePtr,
            DataPathManager::getBigDataDir() / path("pubchem_300000.sdf"));
    }

} // namespace qtr
