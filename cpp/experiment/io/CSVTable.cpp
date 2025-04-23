#include "CSVTable.h"

#include <cassert>

CSVTable::CSVTable(std::vector<std::string> &&columnNames, std::vector<std::vector<std::string>> &&data)
        : _columnNames(std::move(columnNames)), _data(std::move(data)) {
    for (const auto &row: _data) {
        assert(_columnNames.size() == row.size());
    }
}

size_t CSVTable::getColumnCount() const {
    return _columnNames.size();
}

size_t CSVTable::getRowCount() const {
    return _data.size();
}

const std::vector<std::string> &CSVTable::getColumnNames() const {
    return _columnNames;
}

const std::string& CSVTable::getValue(size_t row, size_t col) const {
    return _data.at(row).at(col);
}
