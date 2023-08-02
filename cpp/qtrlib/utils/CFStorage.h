#pragma once

#include <vector>
#include <cstdint>

#include "base_cpp/array.h"
#include "molecule/molecule.h"

namespace qtr {

    class CFStorage {
    public:
        using KeyType = uint64_t;
        using ValueType = indigo::Array<char>;

        ValueType &Add(KeyType key, ValueType &&value);

        CFStorage();

        CFStorage(const CFStorage &) = default;

        CFStorage(CFStorage &&) = default;

        std::unique_ptr<indigo::Molecule> operator[](KeyType key) const;

    private:
        std::vector<ValueType> _map;
    };


} // qtr
