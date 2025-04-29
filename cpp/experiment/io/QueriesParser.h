#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <filesystem>

#include "dataset/SmilesStorage.h"

class QueriesParser {
public:
    explicit QueriesParser(std::filesystem::path  filePath);

    [[nodiscard]] SmilesStorage parse() const;

private:
    std::filesystem::path _filePath;
};
