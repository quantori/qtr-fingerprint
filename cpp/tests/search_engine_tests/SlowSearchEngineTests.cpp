#include "SearchEngineTests.h"
#include "../utils/DataPathManager.h"

#include <gtest/gtest.h>

using namespace qtr;

class SlowTestFixture : public ::testing::Test {
protected:
    std::filesystem::path sdfFile = qtr::DataPathManager::getDataDir() / "pubchem_300000.sdf";
    std::filesystem::path queriesFile = qtr::DataPathManager::getDataDir() / "pubchem_300000_queries.txt";

    void SetUp() override {
        if (DataPathManager::getBigDataDir().empty())
            GTEST_SKIP() << "Skipping all tests with empty big data path";
    }
};

TEST_F(SlowTestFixture, EXHAUSTIVE_BUILD_PUBCHEM_300000) {
    if (!std::filesystem::exists(sdfFile))
        GTEST_SKIP() << "No file " << sdfFile << " to test";
    SearchEngineTests searchEngineTests(SearchEngineFactory::EXHAUSTIVE);
    searchEngineTests.testBuildPubchem300000();
}

TEST_F(SlowTestFixture, EXHAUSTIVE_PUBCHEM_300000) {
    if (!std::filesystem::exists(sdfFile))
        GTEST_SKIP() << "No file " << sdfFile << " to test";
    if (!std::filesystem::exists(queriesFile))
        GTEST_SKIP() << "No file " << queriesFile << " to test";
    SearchEngineTests searchEngineTests(SearchEngineFactory::EXHAUSTIVE);
    searchEngineTests.testPubchem300000();
}

TEST_F(SlowTestFixture, BINGO_BUILD_PUBCHEM_300000) {
    if (!std::filesystem::exists(sdfFile))
        GTEST_SKIP() << "No file " << sdfFile << " to test";
    SearchEngineTests searchEngineTests(SearchEngineFactory::BINGO);
    searchEngineTests.testBuildPubchem300000();
}

TEST_F(SlowTestFixture, BINGO_PUBCHEM_300000) {
    SearchEngineTests searchEngineTests(SearchEngineFactory::BINGO);
    if (!std::filesystem::exists(sdfFile))
        GTEST_SKIP() << "No file " << sdfFile << " to test";
    if (!std::filesystem::exists(queriesFile))
        GTEST_SKIP() << "No file " << queriesFile << " to test";
    searchEngineTests.testPubchem300000();
}

