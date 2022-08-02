#pragma once

#include <iostream>
#include <vector>
#include <filesystem>
#include <fstream>

#include "ColumnsIOConsts.h"

namespace qtr {

    class ColumnsReader {
    public:
        explicit ColumnsReader(std::istream *inStream) : _inStream(inStream) {};

        explicit ColumnsReader(const std::filesystem::path &fileName) : ColumnsReader(new std::ifstream(fileName)) {}

        ColumnsReader(const ColumnsReader&) = delete;

        std::vector<column_t> read();

        ~ColumnsReader();

    private:
        std::istream *_inStream;

    };

} // namespace qtr