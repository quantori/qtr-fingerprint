#pragma once

#include "indigo.h"
#include "IndigoSession.h"
#include "FingerprintTable.h"
#include "IndigoSDFileIterator.h"

template<size_t fingerprintSizeInBytes>
qtr::FingerprintTable<fingerprintSizeInBytes> qtr::buildFingerprintTableFromSDFFile(const std::string &sdfFile) {
    auto indigoSessionPtr = indigo_cpp::IndigoSession::create();
    qtr::FingerprintTable<fingerprintSizeInBytes> table;
    for (auto &mol: indigoSessionPtr->iterateSDFile(sdfFile)) {
        mol->aromatize();
        auto fingerprint = Fingerprint<fingerprintSizeInBytes>(indigoFingerprint(mol->id(), "sub"));
        table.emplace_back(fingerprint);
    }
    return table;
}
