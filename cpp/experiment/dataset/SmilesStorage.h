#pragma once

#include <vector>
#include <string>

class SmilesStorage {
public:

    explicit SmilesStorage(std::vector<std::string> &&smilesSequence) : _smiles(std::move(smilesSequence)) {
    }

    [[nodiscard]] size_t size() const {
        return _smiles.size();
    }

    [[nodiscard]] const std::string &smiles(size_t idx) const {
        return _smiles.at(idx);
    }

private:
    std::vector<std::string> _smiles;
};
