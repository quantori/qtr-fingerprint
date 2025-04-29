#pragma once

#include <unordered_map>
#include "StatRow.h"
#include "CSVTable.h"

class StatTable {
public:
    void addRow(const StatRow &row);

    CSVTable toCSVTable() const;

    const StatRow &operator[](size_t i) const {
        return _rows[i];
    }

    StatRow &operator[](size_t i) {
        return _rows[i];
    }

private:
    std::vector<StatRow> _rows;

    [[nodiscard]] std::unordered_map<std::string, size_t> buildColumnNames() const;

    std::vector<std::vector<std::string>>
    buildCSVData(const std::unordered_map<std::string, size_t> &columnNames) const;
};
