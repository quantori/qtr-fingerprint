#pragma once

#include <filesystem>
#include <vector>
#include <string>
#include <mutex>

class SmilesCSVParser {
public:
    explicit SmilesCSVParser(std::filesystem::path filePath);

    void parse(std::vector<std::string> &dest, std::mutex &m) const;

private:
    std::filesystem::path _filePath;
};

