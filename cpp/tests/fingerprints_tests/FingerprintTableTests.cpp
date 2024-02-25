#include "FingerprintTable.h"

#include "../utils/DataPathManager.h"

#include <gtest/gtest.h>

#include <filesystem>

TEST(FingerprintTable, BUILD) {
    qtr::FingerprintTable table;
}

TEST(FingerprintTable, BuildFromSdf) {
    using namespace std::filesystem;
    std::string sdf = qtr::DataPathManager::getDataDir() / path("pubchem_10.sdf");
    qtr::FingerprintTable table(sdf);
}