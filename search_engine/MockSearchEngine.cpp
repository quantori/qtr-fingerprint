#include "MockSearchEngine.h"

void MockSearchEngine::build(std::string path) {
    std::this_thread::sleep_for(std::chrono::seconds(2));
}

std::vector<indigo_cpp::IndigoMolecule>
MockSearchEngine::findOverMolecules(const indigo_cpp::IndigoQueryMolecule &mol) {
    std::this_thread::sleep_for(std::chrono::milliseconds(15));
    return {};
}

