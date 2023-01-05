#include <algorithm>
#include <cassert>

#include "HuffmanSmilesTable.h"

namespace qtr {

    namespace {
        struct CompareByFirst {
            bool operator()(const std::pair<HuffmanSmilesTable::KeyType, std::vector<bool>> &a,
                            const std::pair<HuffmanSmilesTable::KeyType, std::vector<bool>> &b) const {
                return a.first < b.first;
            }
        };
    }

    HuffmanSmilesTable HuffmanSmilesTable::Builder::build() {
        assert(std::is_sorted(_storage.begin(), _storage.end()));
        return {_storage, _huffmanCoder};
    }

    std::shared_ptr<HuffmanSmilesTable> HuffmanSmilesTable::Builder::buildPtr() {
        assert(std::is_sorted(_storage.begin(), _storage.end()));
        return std::shared_ptr<HuffmanSmilesTable>(new HuffmanSmilesTable(_storage, _huffmanCoder));
    }

    HuffmanSmilesTable::Builder &
    HuffmanSmilesTable::Builder::operator+=(const std::pair<KeyType, std::vector<bool>> &item) {
        _storage.emplace_back(item);
        return *this;
    }

    HuffmanSmilesTable::Builder &HuffmanSmilesTable::Builder::operator+=(const std::pair<KeyType, std::string> &item) {
        return this->operator+=({item.first, _huffmanCoder.encode(item.second)});
    }

    HuffmanSmilesTable::ValueType HuffmanSmilesTable::operator[](const HuffmanSmilesTable::KeyType &key) const {
        auto it = std::lower_bound(_storage.begin(), _storage.end(), std::make_pair(key, std::vector<bool>()),
                                   CompareByFirst());
        assert(it != _storage.end() && it->first == key && "Key wasn't found");
        return _huffmanCoder.decode(it->second);
    }

    HuffmanSmilesTable::ValueType HuffmanSmilesTable::at(const HuffmanSmilesTable::KeyType &key) const {
        auto it = std::lower_bound(_storage.begin(), _storage.end(), std::make_pair(key, std::vector<bool>()),
                                   CompareByFirst());
        if (it == _storage.end() || it->first != key)
            throw std::out_of_range("Key wasn't found in smiles table");
        return _huffmanCoder.decode(it->second);
    }
} // qtr