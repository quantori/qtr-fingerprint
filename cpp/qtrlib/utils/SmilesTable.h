#pragma once

#include <utility>
#include <vector>
#include <cstdint>


#include "HuffmanCoder.h"

namespace qtr {

    class SmilesTable {
    public:
        using KeyType = uint64_t;
        using ValueType = std::vector<bool>;
        using StorageType = std::vector<std::pair<KeyType, ValueType>>;

    private:
        StorageType _storage;
        HuffmanCoder _huffmanCoder;

        SmilesTable(StorageType storage, HuffmanCoder huffmanCoder) : _storage(std::move(storage)),
                                                                      _huffmanCoder(std::move(huffmanCoder)) {}

    public:

        class Builder {
        private:
            StorageType _storage;
            HuffmanCoder _huffmanCoder;
        public:

            Builder(HuffmanCoder huffmanCoder) : _huffmanCoder(std::move(huffmanCoder)) {};

            SmilesTable build();

            Builder &operator+=(const std::pair<KeyType, ValueType> &item);

            Builder &operator+=(const std::pair<KeyType, std::string> &item);


        };

        friend SmilesTable Builder::build();

        SmilesTable() = delete;

        std::string operator[](KeyType key) const;

        [[nodiscard]] std::string at(KeyType key) const;

    };

} // qtr

