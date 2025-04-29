#pragma once

#include <string>
#include <vector>
#include <memory>
#include <fstream>

#include "CSVTable.h"

class CSVWriter {
public:
    CSVWriter();

    explicit CSVWriter(char delimiter);

    bool writeToFile(const std::string &filename, const CSVTable &table, bool includeHeader = true);

    bool writeToStream(std::ostream &stream, const CSVTable &table, bool includeHeader = true);

    void setDelimiter(char delimiter);

private:
    char delimiter_;

    [[nodiscard]] std::string escapeField(const std::string &field) const;
};
