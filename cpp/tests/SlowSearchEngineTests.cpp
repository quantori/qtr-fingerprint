#include "SearchEngineTests.h"
#include "utils/DataPathManager.h"

#include <gtest/gtest.h>

using namespace qtr;

class SlowTestFixture : public ::testing::Test
{
protected:
    void SetUp() override
    {
        ::testing::Test::RecordProperty("SETUP", 1);
        if (DataPathManager::getBigDataDir().empty())
            GTEST_SKIP() << "Skipping all tests with empty big data path";
    }
};

TEST_F(SlowTestFixture, EXHAUSTIVE_PUBCHEM_119697)
{
    SearchEngineTests searchEngineTests(SearchEngineFactory::EXHAUSTIVE);
    searchEngineTests.testPubchem994();
}

// TEST_F(SlowTestFixture, EXHAUSTIVE_BUILD_PUBCHEM_119697)
// {
//     SearchEngineTests searchEngineTests(SearchEngineFactory::EXHAUSTIVE);
//     searchEngineTests.testBuildPubchem119697();
// }


TEST_F(SlowTestFixture, BINGO_PUBCHEM_119697)
{
    SearchEngineTests searchEngineTests(SearchEngineFactory::BINGO);
    searchEngineTests.testPubchem994();
}

// TEST_F(SlowTestFixture, BINGO_BUILD_PUBCHEM_119697)
// {
//     SearchEngineTests searchEngineTests(SearchEngineFactory::BINGO);
//     searchEngineTests.testBuildPubchem119697();
// }