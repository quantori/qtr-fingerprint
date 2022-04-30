#pragma once

#include "SearchEngineInterface.h"

#include <thread>

class MockSearchEngine final : public SearchEngineInterface {
public:
    inline MockSearchEngine() = default;
    inline ~MockSearchEngine() final = default;

    inline void build(std::string path) final {
        std::this_thread::sleep_for(std::chrono::seconds(2));
    }
    inline std::vector<indigo_cpp::IndigoMolecule> findOverMolecules(indigo_cpp::IndigoMolecule mol) override {
        std::this_thread::sleep_for(std::chrono::milliseconds (15));
        return {};
    }
};