#pragma once


#include <cstdint>
#include <cstring>
#include "IndigoMolecule.h"
#include "Utils.h"
#include "indigo.h"

template<std::size_t FINGERPRINT_SIZE>
class Fingerprint;

template<std::size_t FINGERPRINT_SIZE>
Fingerprint<FINGERPRINT_SIZE>::Fingerprint(indigo_cpp::IndigoMolecule &mol) {
    mol.aromatize();
    int fingerprint = indigoFingerprint(mol.id(), "sub");
    const char *fp = indigoToString(fingerprint);
    buildFromIndigoFingerprint(fp);

}

template<std::size_t FINGERPRINT_SIZE>
void Fingerprint<FINGERPRINT_SIZE>::buildFromIndigoFingerprint(const char *fp) {
    for (int i = 0; fp[i] != '\0'; ++i) {
        BlockType intNumber = (fp[i] >= 'a') ? (fp[i] - 'a' + 10) : (fp[i] - '0');
        _data[getBlockNumber(i)] += intNumber << (blockTypeBits - lenInBinary(i + 1));
    }
}