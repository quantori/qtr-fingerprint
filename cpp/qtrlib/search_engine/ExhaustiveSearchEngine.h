#pragma once

#include "FingerprintTableSearchEngine.h"

namespace qtr {

class ExhaustiveSearchEngine : public FingerprintTableSearchEngine {
public:
    using FingerprintTableSearchEngine::FingerprintTableSearchEngine;

    void build(const std::string &path) override;

private:
    std::vector<const IndigoFingerprintTableView *> findTableViews(const qtr::IndigoFingerprint &fp) const override;

    IndigoFingerprintTableView _tableView;
};

} // namespace qtr