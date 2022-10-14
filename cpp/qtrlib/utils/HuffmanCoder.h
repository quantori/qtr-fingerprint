#pragma once

#include <filesystem>
#include <vector>

namespace qtr {
    class HuffmanCoder {
    private:

        struct Node {
            char symbol;
            size_t left, right;

            Node(char symbol) : symbol(symbol), left(-1), right(-1) {}

            Node(char symbol, size_t left, size_t right) : symbol(symbol), left(left), right(right) {}
        };

        std::vector<Node> _treeNodes;
        size_t _treeRoot;
        std::vector<bool> _symbolsMap[std::numeric_limits<uint8_t>::max() + 1];

        void createTree(const std::filesystem::path &prioritiesFile);

        void createSymbolsMap();

        const std::vector<bool> &getSymbolCode(char symbol) const;

        std::vector<bool>& getSymbolCode(char symbol);
    public:

        HuffmanCoder(const std::filesystem::path &prioritiesFile);

        std::vector<bool> encode(const std::string &string) const;

        std::string decode(const std::vector<bool> &code) const;

    };

} // qtr

