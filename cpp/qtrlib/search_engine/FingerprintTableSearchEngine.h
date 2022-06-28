#pragma once

#include "FingerprintTable.h"
#include "FingerprintTableView.h"
#include "SearchEngineInterface.h"

namespace qtr {

class FingerprintTableSearchEngine : public SearchEngineInterface {
public:
    FingerprintTableSearchEngine() = delete;
    explicit FingerprintTableSearchEngine(const indigo_cpp::IndigoSessionPtr &indigoSessionPtr);

    ~FingerprintTableSearchEngine() override;

    void build(const std::string &path) override;

    std::vector<indigo_cpp::IndigoMolecule> findOverMolecules(const indigo_cpp::IndigoQueryMolecule &mol) override;

protected:
    virtual std::vector<const IndigoFingerprintTableView *> findTableViews(const qtr::IndigoFingerprint &fp) const = 0;

    IndigoFingerprintTable _fingerprintTable;
    std::vector<indigo_cpp::IndigoMolecule> _molecules;
    indigo_cpp::IndigoSessionPtr _indigoSessionPtr;
};

} // namespace qtr