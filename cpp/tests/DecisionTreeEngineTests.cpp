#include "Common.h"
#include "SearchEngineTests.h"

#include <gtest/gtest.h>

using namespace qtr;

class DecisionTreeSearchEngineTestFixture : public ::testing::Test {
protected:
    DecisionTreeSearchEngineTestFixture()
        : searchEngineTests(SearchEngineFactory::DECISION_TREE)
    { }

    SearchEngineTests searchEngineTests;
};

TEST_F(DecisionTreeSearchEngineTestFixture, DECISION_TREE_PUBCHEM_10) 
{
    searchEngineTests.testPubchem10();
}

TEST_F(DecisionTreeSearchEngineTestFixture, DECISION_TREE_PUBCHEM_100) 
{
    searchEngineTests.testPubchem100();
}

// Uncomment this test to get memory overflow

//TEST_F(DecisionTreeSearchEngineTestFixture, DECISION_TREE_PUBCHEM_994) 
//{
//    searchEngineTests.testPubchem994();
//}