#pragma once

#include "SmilesTable.h"

#include <map>
#include <vector>

namespace qtr {

    class MapSmilesTable : public SmilesTable {
    private:
        std::map<KeyType, ValueType> _map;
    public:
        template<class InputIterator>
        MapSmilesTable(InputIterator first, InputIterator last) : _map(first, last) {}

        MapSmilesTable() = default;

        ValueType operator[](const KeyType &key) const override;

        ValueType at(const KeyType &key) const override;
    };

} // qtr

