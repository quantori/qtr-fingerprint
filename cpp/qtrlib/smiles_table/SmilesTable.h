#pragma once

#include <cstdint>
#include <string>

namespace qtr {

    class SmilesTable {
    public:
        using KeyType = uint64_t;
        using ValueType = std::string;

        virtual ValueType operator[](const KeyType &key) const = 0;

        [[nodiscard]] virtual ValueType at(const KeyType &key) const = 0;

        virtual ~SmilesTable() = default;
    };

} // qtr

