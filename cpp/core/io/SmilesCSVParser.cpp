#include "SmilesCSVParser.h"

#include <glog/logging.h>
#include <fstream>
#include <mutex>

SmilesCSVParser::SmilesCSVParser(std::filesystem::path filePath) : _filePath(std::move(filePath)) {

}

void SmilesCSVParser::parse(std::vector<std::string> &dest, std::mutex &m) const {
    std::ifstream in(_filePath);
    std::string smiles;
    while (in.peek() != EOF) {
        std::string id;
        in >> id >> smiles;
        std::string lineEnding;
        std::getline(in, lineEnding);
        if (smiles.empty()) {
            LOG(WARNING) << "Found line with wrong formatting and skipped it: \"" << smiles << lineEnding << "\"";
            continue;
        }
        {
            std::lock_guard<std::mutex> lockGuard(m);
            dest.push_back(smiles);
        }
    }
}
