#pragma once

#include <vector>
#include "base_cpp/array.h"

class IndigoQueryFingerprint {
public:
    explicit IndigoQueryFingerprint(const indigo::Array<byte> &fingerprint);

    [[nodiscard]] bool isSubFingerprint(const indigo::Array<byte> &fingerprint) const;

private:
    std::vector<std::pair<int, byte>> _bytes;
};

