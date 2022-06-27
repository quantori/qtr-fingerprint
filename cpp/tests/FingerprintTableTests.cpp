#include "FingerprintTable.h"

#include <gtest/gtest.h>

TEST(FingerprintTable, BUILD) {
    qtr::IndigoFingerprintTable table;
}

TEST(FingerprintTable, BuildFromSdf) {
    qtr::IndigoFingerprintTable table = qtr::buildFingerprintTableFromSDFFile<467>("../../../data/pubchem_10.sdf");
}