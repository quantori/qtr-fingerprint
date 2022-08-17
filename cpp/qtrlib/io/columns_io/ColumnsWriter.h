#pragma once

#include <filesystem>
#include <iostream>
#include <fstream>
#include <vector>

#include "ColumnsIOConsts.h"

namespace qtr {

    class ColumnsWriter {
    public:
        explicit ColumnsWriter(std::ostream *outStream) : _outStream(outStream) {};

        explicit ColumnsWriter(const std::filesystem::path& fileName) : _outStream(new std::ofstream(fileName)) {};

        ColumnsWriter(const ColumnsWriter &) = default;

        void write(const std::vector<size_t> &columns);

        ~ColumnsWriter();

    private:
        std::ostream *_outStream;

    };

} // namespace qtr
