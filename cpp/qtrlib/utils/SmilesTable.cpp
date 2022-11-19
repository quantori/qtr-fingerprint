#include <algorithm>
#include <cassert>

#include "SmilesTable.h"

namespace qtr {

    namespace {
        struct CompareByFirst {
            bool operator()(const std::pair<SmilesTable::KeyType, SmilesTable::ValueType> &a,
                            const std::pair<SmilesTable::KeyType, SmilesTable::ValueType> &b) const {
                return a.first < b.first;
            }
        };
    }

    SmilesTable SmilesTable::Builder::build() {
        std::sort(_storage.begin(), _storage.end(), CompareByFirst());
        return {_storage, _huffmanCoder};
    }

    SmilesTable::Builder &SmilesTable::Builder::operator+=(const std::pair<KeyType, ValueType> &item) {
        _storage.emplace_back(item);
        return *this;
    }

    SmilesTable::Builder &SmilesTable::Builder::operator+=(const std::pair<KeyType, std::string> &item) {
        return this->operator+=({item.first, _huffmanCoder.encode(item.second)});
    }

    std::string SmilesTable::operator[](SmilesTable::KeyType key) const {
        auto it = std::lower_bound(_storage.begin(), _storage.end(), std::make_pair(key, ValueType()),
                                   CompareByFirst());
        assert(it != _storage.end() && it->first == key && "Key wasn't found");
        return _huffmanCoder.decode(it->second);
    }

    std::string SmilesTable::at(SmilesTable::KeyType key) const {
        auto it = std::lower_bound(_storage.begin(), _storage.end(), std::make_pair(key, ValueType()),
                                   CompareByFirst());
        if (it == _storage.end() || it->first != key)
            throw std::out_of_range("Key wasn't found in smiles table");
        return _huffmanCoder.decode(it->second);
    }
} // qtr