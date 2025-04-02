#include "QueriesParser.h"

#include <utility>

#include "glog/logging.h"

QueriesParser::QueriesParser(std::filesystem::path filePath) : _filePath(std::move(filePath)) {
}

SmilesStorage QueriesParser::parse() const {
    std::ifstream in(_filePath);
    std::vector<std::string> result;
    while (in.peek() != EOF) {
        std::string smiles;
        in >> smiles;
        std::string lineEnding;
        std::getline(in, lineEnding);
        if (smiles.empty()) {
            LOG(ERROR) << "Wrong formatting in queries file";
        }
        result.push_back(smiles);
    }
    LOG(INFO) << "Loaded " << result.size() << " molecules from file " << _filePath;
    return SmilesStorage(std::move(result));
}
