#pragma once

#include <vector>

#include "StatEntry.h"

struct StatRow {
    template<typename T>
    void addEntry(const std::string &name, const T &value) {
        _entries.emplace_back(name, value);
    }

    inline StatEntry &operator[](size_t i) {
        return _entries[i];
    }

    inline const StatEntry &operator[](size_t i) const {
        return _entries[i];
    }

    [[nodiscard]] inline size_t size() const {
        return _entries.size();
    }

private:
    std::vector<StatEntry> _entries;
};
