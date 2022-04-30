#include "MockSearchEngine.h"

#include <thread>

void MockSearchEngine::build(std::string path) {
    std::this_thread::sleep_for(std::chrono::seconds(2));
}

std::vector<indigo_cpp::IndigoMolecule> MockSearchEngine::findOverMolecules(indigo_cpp::IndigoMolecule mol) {
    std::this_thread::sleep_for(std::chrono::milliseconds (15));
    return {};
}

