#pragma once

#include <utility>
#include <vector>
#include <cstdint>

#include "SmilesTable.h"
#include "HuffmanCoder.h"

namespace qtr {

    class HuffmanSmilesTable : public SmilesTable {
    private:
        using StorageType = std::vector<std::pair<KeyType, std::vector<bool>>>;

        StorageType _storage;
        HuffmanCoder _huffmanCoder;

        HuffmanSmilesTable(StorageType storage, HuffmanCoder huffmanCoder) : _storage(std::move(storage)),
                                                                             _huffmanCoder(std::move(huffmanCoder)) {}

    public:
        class Builder {
        private:
            StorageType _storage;
            HuffmanCoder _huffmanCoder;
        public:

            explicit Builder(HuffmanCoder huffmanCoder) : _huffmanCoder(std::move(huffmanCoder)) {};

            HuffmanSmilesTable build();

            std::shared_ptr<HuffmanSmilesTable> buildPtr();

            Builder &operator+=(const std::pair<KeyType, std::vector<bool>> &item);

            Builder &operator+=(const std::pair<KeyType, std::string> &item);


        };

        friend HuffmanSmilesTable Builder::build();

        HuffmanSmilesTable() = delete;

        ValueType operator[](const KeyType &key) const override;

        [[nodiscard]] ValueType at(const KeyType &key) const override;

    };

} // qtr

