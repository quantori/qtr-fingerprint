#include "SmilesDirParser.h"

#include "SmilesCSVParser.h"

#include <glog/logging.h>

SmilesDirParser::SmilesDirParser(std::filesystem::path dirPath) : _dirPath(std::move(dirPath)) {
}

std::vector<std::string> SmilesDirParser::parse() const {
    std::vector<std::string> results;
    for (auto &dirEntry: std::filesystem::directory_iterator(_dirPath)) {
        auto &filename = dirEntry.path();
        if (!std::filesystem::is_regular_file(filename)) {
            LOG(WARNING) << "Skipped entity: " << filename;
            continue;
        }
        auto fileResults = SmilesCSVParser(filename).parse();
        results.insert(results.end(), fileResults.begin(), fileResults.end());
    }
    return results;
}


