#include "common.h"
#include "SearchEngineFactory.h"

#include <gtest/gtest.h>

#include <filesystem>

using namespace indigo_cpp;
using namespace std::filesystem;


class BingoSearchEngineTestFixture : public ::testing::Test {
protected:
    BingoSearchEngineTestFixture()
        : _indigoSessionPtr(IndigoSession::create())
        , _searchEnginePtr(SearchEngineFactory::create(_indigoSessionPtr))
        , _dataDir(getDataDir())
    { }

    IndigoSessionPtr _indigoSessionPtr;
    SearchEnginePtr _searchEnginePtr;
    path _dataDir;
};

TEST_F(BingoSearchEngineTestFixture, BINGO_PUBCHEM_10) 
{
    testSearchEngine(
        _searchEnginePtr,
        _indigoSessionPtr,
        _dataDir / path("pubchem_10.sdf"),
        _dataDir / path("pubchem_10_queries.txt")
    );
}

TEST_F(BingoSearchEngineTestFixture, BINGO_PUBCHEM_100) 
{
    testSearchEngine(
        _searchEnginePtr,
        _indigoSessionPtr,
        _dataDir / path("pubchem_100.sdf"),
        _dataDir / path("pubchem_100_queries.txt")
    );
}

TEST_F(BingoSearchEngineTestFixture, BINGO_PUBCHEM_994) 
{
    testSearchEngine(
        _searchEnginePtr,
        _indigoSessionPtr,
        _dataDir / path("pubchem_994.sdf"),
        _dataDir / path("pubchem_994_queries.txt")
    );
}