#pragma once

#include "FingerprintTable.h"
#include "TableView.h"

#include <memory>
#include <vector>

namespace qtr {

template<size_t fingerprintSizeInBytes>
using FingerprintTableView = TableView<FingerprintTable<fingerprintSizeInBytes>>;

using IndigoFingerprintTableView = FingerprintTableView<IndigoFingerprint::sizeInBytes>;

} // namespace qtr