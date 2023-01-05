#include "MapSmilesTable.h"

#include <stdexcept>

namespace qtr {
    SmilesTable::ValueType MapSmilesTable::operator[](const SmilesTable::KeyType &key) const {
        return _map.find(key)->second;
    }

    SmilesTable::ValueType MapSmilesTable::at(const SmilesTable::KeyType &key) const {
        auto it = _map.find(key);
        if (it == _map.end())
            throw std::out_of_range("Key wasn't found in smiles table");
        return it->second;
    }

    MapSmilesTable::MapSmilesTable(const std::vector<std::pair<KeyType, ValueType>> &values)
            : _map(values.begin(), values.end()) {}

} // qtr