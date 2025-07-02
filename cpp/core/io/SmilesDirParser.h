#pragma once

#include <filesystem>

#include "dataset/SmilesStorage.h"

class SmilesDirParser {
public:
    explicit SmilesDirParser(std::filesystem::path dirPath);

    [[nodiscard]] SmilesStorage parse() const;
private:
    std::filesystem::path _dirPath;
};
