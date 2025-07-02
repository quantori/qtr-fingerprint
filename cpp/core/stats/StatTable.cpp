#include "StatTable.h"

#include <unordered_set>

void StatTable::addRow(const StatRow &row) {
    _rows.emplace_back(row);
}

CSVTable StatTable::toCSVTable() const {
    auto numberedColumnNames = buildColumnNames();
    auto data = buildCSVData(numberedColumnNames);
    std::vector<std::string> columnNames(numberedColumnNames.size());
    for (auto &[columnName, idx]: numberedColumnNames) {
        columnNames[idx] = columnName;
    }
    return {std::move(columnNames), std::move(data)};
}

std::unordered_map<std::string, size_t> StatTable::buildColumnNames() const {
    std::unordered_map<std::string, size_t> columnNames;
    for (auto &row: _rows) {
        for (size_t i = 0; i < row.size(); i++) {
            auto &name = row[i].name;
            if (columnNames.contains(name)) {
                continue;
            }
            columnNames[name] = columnNames.size();
        }
    }
    return columnNames;
}

std::vector<std::vector<std::string>>
StatTable::buildCSVData(const std::unordered_map<std::string, size_t> &columnNames) const {
    std::vector<std::vector<std::string>> data;
    for (auto &row: _rows) {
        std::vector<std::string> dataRow(columnNames.size(), "");
        for (size_t i = 0; i < row.size(); i++) {
            size_t columnIdx = columnNames.at(row[i].name);
            dataRow[columnIdx] = row[i].value;
        }
        data.emplace_back(dataRow);
    }
    return data;
}
