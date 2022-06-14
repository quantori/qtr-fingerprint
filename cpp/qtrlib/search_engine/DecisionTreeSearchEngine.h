#pragma once

#include "FingerprintTable.h"
#include "SearchEngineInterface.h"

#include "IndigoMolecule.h"

#include <vector>

namespace qtr {

class DecisionTreeSearchEngine : public SearchEngineInterface {
public:
    DecisionTreeSearchEngine() = delete;
    explicit DecisionTreeSearchEngine(const indigo_cpp::IndigoSessionPtr &indigoSessionPtr);

    ~DecisionTreeSearchEngine() override;

    void build(const std::string &path) override;

    std::vector<indigo_cpp::IndigoMolecule> findOverMolecules(const indigo_cpp::IndigoQueryMolecule &mol) override;

private:
    IndigoFingerprintTable _fingerprintTable;
    std::vector<indigo_cpp::IndigoMolecule> _molecules;
    indigo_cpp::IndigoSessionPtr _indigoSessionPtr;
};

} // namespace qtr