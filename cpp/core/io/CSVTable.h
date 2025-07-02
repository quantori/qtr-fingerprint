#pragma once

#include <vector>
#include <string>

class CSVTable {
public:
    CSVTable(std::vector<std::string> &&columnNames, std::vector<std::vector<std::string>> &&data);

    [[nodiscard]] size_t getColumnCount() const;

    [[nodiscard]] size_t getRowCount() const;

    [[nodiscard]] const std::vector<std::string> &getColumnNames() const;

    [[nodiscard]] const std::string &getValue(size_t row, size_t col) const;

private:
    std::vector<std::string> _columnNames;
    std::vector<std::vector<std::string>> _data;
};