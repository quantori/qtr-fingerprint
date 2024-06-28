#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <filesystem>

class QueriesParser {
public:
    explicit QueriesParser(std::filesystem::path  filePath);

    [[nodiscard]] std::vector<std::string> parse() const;

private:
    std::filesystem::path _filePath;
};
