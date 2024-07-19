#pragma once

#include <filesystem>
#include <string>
#include <vector>
#include <mutex>

class SmilesDirParser {
public:
    explicit SmilesDirParser(std::filesystem::path dirPath);

    void parse(std::vector<std::string> &dest, std::mutex &m) const;
private:
    std::filesystem::path _dirPath;
};
