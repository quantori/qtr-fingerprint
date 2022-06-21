#include "ExhaustiveSearchEngine.h"
#include "QtrIndigoFingerprint.h"
#include "Utils.h"

#include "IndigoSDFileIterator.h"
#include "IndigoSubstructureMatcher.h"
#include "indigo.h"

#include <glog/logging.h>

using namespace indigo_cpp;

namespace qtr {

void ExhaustiveSearchEngine::build(const std::string &path)
{
    FingerprintTableSearchEngine::build(path);
    _tableView = IndigoFingerprintTableView(&_fingerprintTable);
}

std::vector<const IndigoFingerprintTableView *> ExhaustiveSearchEngine::findTableViews(const qtr::IndigoFingerprint &fp) const
{
    return {&_tableView};
}

} // namespace qtr