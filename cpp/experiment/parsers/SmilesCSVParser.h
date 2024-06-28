#pragma once

#include <filesystem>
#include <vector>
#include <string>

class SmilesCSVParser {
public:
    explicit SmilesCSVParser(std::filesystem::path filePath);

    std::vector<std::string> parse() const;

private:
    std::filesystem::path _filePath;
};

