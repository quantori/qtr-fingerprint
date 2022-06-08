#pragma once

#include "SearchEngineInterface.h"

#include "QtrIndigoFingerprint.h"
#include "IndigoMolecule.h"

class ExhaustiveSearchEngine : public SearchEngineInterface {
public:
    ExhaustiveSearchEngine() = delete;
    explicit ExhaustiveSearchEngine(const indigo_cpp::IndigoSessionPtr &indigoSessionPtr);

    ~ExhaustiveSearchEngine() override;

    void build(const std::string &path) override;

    std::vector<indigo_cpp::IndigoMolecule> findOverMolecules(const indigo_cpp::IndigoQueryMolecule &mol) override;

private:

    struct MoleculeInfo {
        indigo_cpp::IndigoMolecule molecule;
        QtrIndigoFingerprint fingerprint;
    };

    std::vector<MoleculeInfo> _moleculesInfo;
    indigo_cpp::IndigoSessionPtr _indigoSessionPtr;
};
