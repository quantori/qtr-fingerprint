#include "Fingerprint.h"

#include "QtrIndigoFingerprint.h"

namespace qtr {

template<>
IndigoFingerprint::Fingerprint(const QtrIndigoFingerprint &f) {
    setBytes(f.data());
}

} // namespace qtr