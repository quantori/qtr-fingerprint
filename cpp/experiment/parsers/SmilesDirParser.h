#pragma once

#include <filesystem>
#include <string>
#include <vector>

class SmilesDirParser {
public:
    explicit SmilesDirParser(std::filesystem::path dirPath);

    std::vector<std::string> parse() const;
private:
    std::filesystem::path _dirPath;
};
