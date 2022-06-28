#include "Common.h"
#include "SearchEngineTests.h"

#include <gtest/gtest.h>

using namespace qtr;

class DecisionTreeTrivialSearchEngineTestFixture : public ::testing::Test {
protected:
    DecisionTreeTrivialSearchEngineTestFixture()
        : searchEngineTests(SearchEngineFactory::DECISION_TREE_TRIVIAL)
    { }

    SearchEngineTests searchEngineTests;
};

TEST_F(DecisionTreeTrivialSearchEngineTestFixture, DECISION_TREE_TRIVIAL_PUBCHEM_10) 
{
    searchEngineTests.testPubchem10();
}

TEST_F(DecisionTreeTrivialSearchEngineTestFixture, DECISION_TREE_TRIVIAL_PUBCHEM_100) 
{
    searchEngineTests.testPubchem100();
}

TEST_F(DecisionTreeTrivialSearchEngineTestFixture, DECISION_TREE_TRIVIAL_PUBCHEM_994) 
{
   searchEngineTests.testPubchem994();
}

/////////////////////////////////////////////////////////////////////////////////////

class DecisionTreeOptimalSearchEngineTestFixture : public ::testing::Test {
protected:
    DecisionTreeOptimalSearchEngineTestFixture()
        : searchEngineTests(SearchEngineFactory::DECISION_TREE_OPTIMAL)
    { }

    SearchEngineTests searchEngineTests;
};

TEST_F(DecisionTreeOptimalSearchEngineTestFixture, DECISION_TREE_OPTIMAL_PUBCHEM_10) 
{
    searchEngineTests.testPubchem10();
}

TEST_F(DecisionTreeOptimalSearchEngineTestFixture, DECISION_TREE_OPTIMAL_PUBCHEM_100) 
{
    searchEngineTests.testPubchem100();
}

TEST_F(DecisionTreeOptimalSearchEngineTestFixture, DECISION_TREE_OPTIMAL_PUBCHEM_994) 
{
   searchEngineTests.testPubchem994();
}
