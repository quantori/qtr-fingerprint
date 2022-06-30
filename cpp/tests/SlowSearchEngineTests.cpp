#include "SearchEngineTests.h"
#include "utils/DataPathManager.h"

#include <gtest/gtest.h>

using namespace qtr;

class SlowTestFixture : public ::testing::Test
{
protected:
    void SetUp() override
    {
        if (DataPathManager::getBigDataDir().empty())
            GTEST_SKIP() << "Skipping all tests with empty big data path";
    }
};

TEST_F(SlowTestFixture, EXHAUSTIVE_BUILD_PUBCHEM_300000)
{
    SearchEngineTests searchEngineTests(SearchEngineFactory::EXHAUSTIVE);
    searchEngineTests.testBuildPubchem300000();
}

TEST_F(SlowTestFixture, EXHAUSTIVE_PUBCHEM_300000)
{
    SearchEngineTests searchEngineTests(SearchEngineFactory::EXHAUSTIVE);
    searchEngineTests.testPubchem300000();
}

TEST_F(SlowTestFixture, BINGO_BUILD_PUBCHEM_300000)
{
    SearchEngineTests searchEngineTests(SearchEngineFactory::BINGO);
    searchEngineTests.testBuildPubchem300000();
}

TEST_F(SlowTestFixture, BINGO_PUBCHEM_300000)
{
    SearchEngineTests searchEngineTests(SearchEngineFactory::BINGO);
    searchEngineTests.testPubchem300000();
}

TEST_F(SlowTestFixture, DECISION_TREE_TRIVIAL_BUILD_PUBCHEM_300000)
{
    SearchEngineTests searchEngineTests(SearchEngineFactory::DECISION_TREE_TRIVIAL);
    searchEngineTests.testBuildPubchem300000();
}

TEST_F(SlowTestFixture, DECISION_TREE_TRIVIAL_PUBCHEM_300000)
{
    SearchEngineTests searchEngineTests(SearchEngineFactory::DECISION_TREE_TRIVIAL);
    searchEngineTests.testPubchem300000();
}

TEST_F(SlowTestFixture, DECISION_TREE_OPTIMAL_BUILD_PUBCHEM_300000)
{
    SearchEngineTests searchEngineTests(SearchEngineFactory::DECISION_TREE_OPTIMAL);
    searchEngineTests.testBuildPubchem300000();
}

TEST_F(SlowTestFixture, DECISION_TREE_OPTIMAL_PUBCHEM_300000)
{
    SearchEngineTests searchEngineTests(SearchEngineFactory::DECISION_TREE_OPTIMAL);
    searchEngineTests.testPubchem300000();
}
