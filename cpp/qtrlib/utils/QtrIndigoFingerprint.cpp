#include "QtrIndigoFingerprint.h"

#include "indigo.h"

#include <cstring>

using namespace indigo_cpp;

QtrIndigoFingerprint::QtrIndigoFingerprint(const IndigoBaseMolecule &molecule, const std::string &type)
        : IndigoObject(indigoFingerprint(molecule.id(), type.c_str()), molecule.session()) {}

std::vector<std::byte> QtrIndigoFingerprint::data() const {
    char *data = nullptr;
    int size = 0;

    session()->_checkResult(indigoToBuffer(id(), &data, &size));

    const int zeroesBeg = std::min(size, 203);
    const int zeroesEnd = std::min(size, 347);

    std::vector<std::byte> result(size - zeroesEnd + zeroesBeg);
    memcpy(result.data(), data, zeroesBeg);
    memcpy(result.data() + zeroesBeg, data + zeroesEnd, size - zeroesEnd);

    return result;
}
