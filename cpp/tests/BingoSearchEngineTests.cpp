#include "Common.h"
#include "SearchEngineTests.h"

#include <gtest/gtest.h>

using namespace qtr;

class BingoSearchEngineTestFixture : public ::testing::Test {
protected:
    BingoSearchEngineTestFixture()
        : searchEngineTests(SearchEngineFactory::BINGO)
    { }

    SearchEngineTests searchEngineTests;
};

TEST_F(BingoSearchEngineTestFixture, BINGO_PUBCHEM_10) 
{
    searchEngineTests.testPubchem10();
}

TEST_F(BingoSearchEngineTestFixture, BINGO_PUBCHEM_100) 
{
    searchEngineTests.testPubchem100();
}

TEST_F(BingoSearchEngineTestFixture, BINGO_PUBCHEM_994) 
{
    searchEngineTests.testPubchem994();
}