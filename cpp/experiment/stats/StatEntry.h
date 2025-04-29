#pragma once

#include <string>
#include <sstream>

struct StatEntry {
    std::string name;
    std::string value;

    template<typename T>
    StatEntry(std::string name, const T &val): name(std::move(name)) {
        std::ostringstream oss;
        oss << val;
        value = oss.str();
    }
};