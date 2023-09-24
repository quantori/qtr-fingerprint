#pragma once

#include <vector>
#include <map>
#include <string>
#include <numeric>
#include <filesystem>

namespace qtr {
    class HuffmanCoder {
    private:

        struct Node {
            char symbol;
            size_t left, right;

            Node(char symbol) : symbol(symbol), left(-1), right(-1) {}

            Node(char symbol, size_t left, size_t right) : symbol(symbol), left(left), right(right) {}
        };

        std::map<char, uint64_t> _symbolsFrequency;
        std::vector<Node> _treeNodes;
        size_t _treeRoot;
        std::vector<bool> _symbolsMap[std::numeric_limits<uint8_t>::max() + 1];

        void createTree(const std::map<char, uint64_t> &symbolsFrequency);

        void createSymbolsMap();

        const std::vector<bool> &getSymbolCode(char symbol) const;

        std::vector<bool> &getSymbolCode(char symbol);

        HuffmanCoder(const std::map<char, uint64_t> &symbolsFrequency);

    public:

        class Builder {
        private:
            std::map<char, uint64_t> symbolsFrequency;

        public:
            Builder() = default;

            Builder &operator+=(const std::string &s);

            Builder &operator+=(const Builder &other);

            HuffmanCoder build() const;
        };

        friend HuffmanCoder Builder::build() const;

        HuffmanCoder() = delete;

        std::vector<bool> encode(const std::string &string) const;

        std::string decode(const std::vector<bool> &code) const;

        void dump(const std::filesystem::path &filePath) const;

        static HuffmanCoder load(const std::filesystem::path &filePath);

    };

} // qtr

