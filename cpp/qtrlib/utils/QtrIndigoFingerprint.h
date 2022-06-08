#pragma once

#include "IndigoBaseMolecule.h"

#include <memory>
#include <string>

class QtrIndigoFingerprint : public indigo_cpp::IndigoObject {
public:
    QtrIndigoFingerprint(const indigo_cpp::IndigoBaseMolecule &molecule, const std::string &type);
    ~QtrIndigoFingerprint() final = default;

    QtrIndigoFingerprint(const QtrIndigoFingerprint &) = delete;
    QtrIndigoFingerprint & operator=(const QtrIndigoFingerprint &) = delete;

    QtrIndigoFingerprint(QtrIndigoFingerprint&&) = default;
    QtrIndigoFingerprint& operator=(QtrIndigoFingerprint&&) = default;

    int countBits() const;
    static int commonBits(const QtrIndigoFingerprint &f1, const QtrIndigoFingerprint &f2);
};

using QtrIndigoFingerprintSPtr = std::shared_ptr<QtrIndigoFingerprint>;