#pragma once

#include "SmilesTable.h"

#include <map>
#include <vector>

namespace qtr {

    class MapSmilesTable : public SmilesTable {
    private:
        std::map<KeyType, ValueType> _map;
    public:
        MapSmilesTable(const std::vector<std::pair<KeyType, ValueType>>& values);

        ValueType operator[](const KeyType &key) const override;

        ValueType at(const KeyType &key) const override;
    };

} // qtr

