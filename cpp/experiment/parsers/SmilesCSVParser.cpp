#include "SmilesCSVParser.h"

#include <glog/logging.h>
#include <fstream>

SmilesCSVParser::SmilesCSVParser(std::filesystem::path filePath) : _filePath(std::move(filePath)) {

}

std::vector<std::string> SmilesCSVParser::parse() const {
    std::ifstream in(_filePath);
    std::string smiles;
    std::vector<std::string> results;
    while (in.peek() != EOF) {
        std::string id;
        in >> id >> smiles;
        std::string lineEnding;
        std::getline(in, lineEnding);
        if (smiles.empty()) {
            LOG(WARNING) << "Found line with wrong formatting and skipped it: \"" << smiles << lineEnding << "\"";
            continue;
        }
        results.push_back(smiles);
    }
    return results;
}
