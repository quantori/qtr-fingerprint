#include "QtrIndigoFingerprint.h"

#include "indigo.h"

using namespace indigo_cpp;

QtrIndigoFingerprint::QtrIndigoFingerprint(const IndigoBaseMolecule &molecule, const std::string &type)
    : IndigoObject(indigoFingerprint(molecule.id(), type.c_str()), molecule.session())
{}

int QtrIndigoFingerprint::countBits() const
{
    return session()->_checkResult(indigoCountBits(id()));
}

int QtrIndigoFingerprint::commonBits(const QtrIndigoFingerprint &f1, const QtrIndigoFingerprint &f2)
{
    return f1.session()->_checkResult(indigoCommonBits(f1.id(), f2.id()));
}