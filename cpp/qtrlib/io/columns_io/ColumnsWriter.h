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

        explicit ColumnsWriter(std::filesystem::path fileName);

        ColumnsWriter(const ColumnsWriter &) = default;

        void write(const std::vector<column_t> &columns);

        ~ColumnsWriter();

    private:
        std::ostream *_outStream;

    };

} // namespace qtr