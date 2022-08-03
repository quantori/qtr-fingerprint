#include "../Common.h"
#include "SearchEngineTests.h"
#include <gtest/gtest.h>

using namespace qtr;

class ExhastiveSearchEngineTestFixture : public ::testing::Test
{
protected:
    ExhastiveSearchEngineTestFixture()
        : searchEngineTests(SearchEngineFactory::EXHAUSTIVE)
    {
    }

    SearchEngineTests searchEngineTests;
};

TEST_F(ExhastiveSearchEngineTestFixture, EXHAUSTIVE_PUBCHEM_10)
{
    searchEngineTests.testPubchem10();
}

TEST_F(ExhastiveSearchEngineTestFixture, EXHAUSTIVE_PUBCHEM_100)
{
    searchEngineTests.testPubchem100();
}

TEST_F(ExhastiveSearchEngineTestFixture, EXHAUSTIVE_PUBCHEM_994)
{
    searchEngineTests.testPubchem994();
}