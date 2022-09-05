#pragma once

#include "FingerprintTable.h"
#include "TableView.h"

#include <memory>
#include <vector>

namespace qtr {

template<size_t fingerprintSizeInBytes>
using FingerprintTableView = TableView<FingerprintTable<IndigoFingerprint::size>>;

using IndigoFingerprintTableView = FingerprintTableView<IndigoFingerprint::size>;

} // namespace qtr