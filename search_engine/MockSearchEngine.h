#pragma once

#include "SearchEngineInterface.h"

class MockSearchEngine final : public SearchEngineInterface {
public:
    inline MockSearchEngine() = default;
    inline ~MockSearchEngine() final = default;

    void build(std::string path) override;
    std::vector<indigo_cpp::IndigoMolecule> findOverMolecules(indigo_cpp::IndigoMolecule mol) override;
};