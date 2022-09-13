#include "../Common.h"

#include "Fingerprint.h"
#include "QtrIndigoFingerprint.h"

#include "IndigoMolecule.h"
#include "IndigoSession.h"

#include <gtest/gtest.h>

#include <memory>

using namespace indigo_cpp;
using namespace qtr;

class FingerprintTestFixture : public ::testing::Test {
protected:
    FingerprintTestFixture() {
        IndigoSessionPtr indigoSessionPtr = IndigoSession::create();
        IndigoMolecule molecule = indigoSessionPtr->loadMolecule("O=C(C)Oc1ccccc1C(=O)O");
        QtrIndigoFingerprint fingerprint(molecule, "sub");
        _data = fingerprint.data();
    }

    std::vector<std::byte> _data;
};

TEST_F(FingerprintTestFixture, SIZE) {
    EXPECT_EQ(_data.size(), qtr::IndigoFingerprint::size());
}

TEST_F(FingerprintTestFixture, CONVERSION) {
    auto fingerprint = std::make_unique<qtr::IndigoFingerprint>();
    fingerprint->setBytes(_data);
    std::vector<std::byte> data = fingerprint->getBytes();
    compareTwoVectors(data, _data);
}
