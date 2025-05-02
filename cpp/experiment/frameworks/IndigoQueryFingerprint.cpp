#include "IndigoQueryFingerprint.h"


bool IndigoQueryFingerprint::isSubFingerprint(const indigo::Array<byte> &fingerprint) const {
    for (auto& [idx, b] : _bytes) {
        byte val = fingerprint.at(idx);
        if ((b & val) != b) {
            return false;
        }
    }
    return true;
}

IndigoQueryFingerprint::IndigoQueryFingerprint(const indigo::Array<byte> &fingerprint) {
    for (int i = 0; i < fingerprint.size(); i++) {
        byte b = fingerprint.at(i);
        if (b == 0) {
            continue;
        }
        _bytes.emplace_back(i, b);
    }
}
